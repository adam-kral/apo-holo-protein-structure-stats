import itertools
from typing import Iterator, List

from .base_analyses import CachedAnalyzer, Analyzer

from Bio.PDB import is_aa, NeighborSearch
from Bio.PDB import Entity, Residue, Chain, Model


def get_hetero_atom_residues(struct: Entity) -> Iterator[Residue]:
    """ non-polymer ligands, excl. water

    BioPython's residue is a group of atoms belonging to the same molecule (made possible by mmcifs author_seq_id of heteroatoms,
    label_seq_id is '.' because a heterocompound does not belong to the poly-chain)
    scans all chains (we deliberately ignore author chain assignment in case of a heteroatom)"""

    return filter(lambda res: res.id[0].startswith('H_'), struct.get_residues())  # .id[0] is a residue's hetatm flag


def get_short_peptide_ligands(struct: Entity, peptide_length_limit: int) -> Iterator[Chain]:
    return filter(lambda chain: sum(is_aa(residue) for residue in chain) <= peptide_length_limit, struct.get_chains())


def get_all_ligands(struct: Model) -> Iterator[Entity]:
    return itertools.chain(get_hetero_atom_residues(struct), get_short_peptide_ligands(struct, 15))


class GetChains(CachedAnalyzer):
    def run(self, struct: Model) -> List[Chain]:
        return list(filter(lambda chain: sum(is_aa(residue) for residue in chain) >= 50, struct.get_chains()))


class GetMainChain(CachedAnalyzer):
    def run(self, struct: Model, get_chains: GetChains) -> Chain:
        """ :returns one chain (the first of GetChains), because hopefully we're working with single-chain structures """
        return get_chains(struct)[0]


class IsHolo(CachedAnalyzer):
    def run(self, struct: Model, get_main_chain: GetMainChain):
        def has_at_least_n_non_hydrogen_atoms(ligand, n):
            non_hydrogen_atoms = 0

            for atom in ligand.get_atoms():
                assert atom.element is not None
                if atom.element != 'H':
                    non_hydrogen_atoms += 1

                if non_hydrogen_atoms >= n:
                    return True  # todo nakonec můžu asi sumovat všechny, stejně budu chtít konfigurovatelny output, aby mi dal počet atomů ligandu, nebo budu dělat statistiky, kolik atomu má průměrný ligand atp.

            return False

        # ligand has >= 6 non-hydrogen atoms
        ligands = list(filter(lambda lig: has_at_least_n_non_hydrogen_atoms(lig, 6), get_all_ligands(struct)))

        # ligand is within RADIUS in contact with MIN_RESIDUES_WITHIN_LIGAND residues

        # (in original paper they used a program LPC, for ensuring specific interaction of ligand with at least 6 residue, this is a "shortcut",
        #    a temporary condition (simple))

        main_chain_atoms = list(get_main_chain(struct).get_atoms())
        ns = NeighborSearch(main_chain_atoms)

        RADIUS = 4.5
        MIN_RESIDUES_WITHIN_LIGAND = 4

        acceptable_ligands = []

        for ligand in ligands:
            residues_in_contact_with_ligand = set()
            for ligand_atom in ligand.get_atoms():  # ligand can be a chain or a residue

                main_chain_atoms_in_contact = ns.search(ligand_atom.get_coord(), RADIUS)

                for atom in main_chain_atoms_in_contact:
                    residues_in_contact_with_ligand.add(atom.get_parent())

            if len(residues_in_contact_with_ligand) >= MIN_RESIDUES_WITHIN_LIGAND:
                acceptable_ligands.append(ligand)

        return len(acceptable_ligands) > 0


