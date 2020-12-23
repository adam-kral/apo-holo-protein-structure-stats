import itertools
from typing import Iterator, List

import numpy as np
import rmsd
from Bio.PDB import is_aa, NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Residue import Residue

from .base_analyses import CachedAnalyzer, Analyzer, SerializableCachedAnalyzer


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


def chain_to_polypeptide(chain):
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(chain, aa_only=0)  # allow non-standard aas?

    if len(polypeptides) != 1:
        print('warning ', len(polypeptides), ' polypeptides from one chain, extending first pp')

        for pp in polypeptides[1:]:
            polypeptides[0].extend(pp)

    return polypeptides[0]


# debug aligner for preliminary sequence analysis (apo-holo/holo-holo), to see how they differ, if they differ
from Bio import Align
aligner = Align.PairwiseAligner(mode='global',
                                open_gap_score=-0.5,
                                extend_gap_score=-0.1,
                                end_gap_score=0,)  # match by default 1 and mismatch 0


def get_c_alpha_coords(polypeptide):
    return np.array([res['CA'].get_coord() for res in polypeptide])



class RMSD(SerializableCachedAnalyzer):
    """nebo GetRMSD - naming. todo bude brát nejspíš numpy array, nebo nějakej jeho wrapper"""

    def run(self, struct1: Entity, struct2: Entity, get_main_chain) -> float:
        pp1, pp2 = map(lambda struct: chain_to_polypeptide(get_main_chain(struct)), (struct1, struct2))
        seq1, seq2 = map(lambda pp: pp.get_sequence(), (pp1, pp2))

        if seq1 != seq2:
            # debug print to see how seqs differ

            #     residues_mapping = # might be a fallback to the pdb-residue-to-uniprot-residue mapping api, depends, what we want
            #           (also we have some mapping (segment to segment) when we do get_best_isoform, there can probably be mutations though)
            #
            #          if there are some extra residues on ends, we can truncate the sequences
            #          maybe I wouldn't use the API, just do this alignment, if we wanted to compare structs with non-100%-identical sequences)

            alignment = next(aligner.align(seq1, seq2))
            print(alignment)
            import math
            return math.inf

        P, Q = map(get_c_alpha_coords, (pp1, pp2))

        return rmsd.kabsch_rmsd(P, Q)


