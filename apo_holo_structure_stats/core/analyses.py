from bisect import bisect_left
from enum import Enum, auto
import itertools
from dataclasses import dataclass
from typing import Iterator, List, Dict, Tuple, Any

import freesasa as freesasa
import numpy as np
import rmsd
from Bio.PDB import is_aa, NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Residue import Residue

from apo_holo_structure_stats.input.download import get_secondary_structure
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







class SSType(Enum):
    HELIX = auto()
    STRAND = auto()
    NO_TYPE = auto()  # no ss type found for a residue


@dataclass
class SSForChain:
    """ Can determine if a residue is a secondary structure by bisecting the ss_start lists.

    <ss>_start and <ss>_end lists have to be sorted and for a ss type their positions correspond (e.g. helices_start[i] and helices_end[i]
    mark the start and end of i-th helical segment)

    If we were always to compare the whole sequences of residues, it would be faster to implement it with a pointer to
    current SS start/end segment. But this should be fast enough, is easier to implement (just using bisect module),
    and more versatile use (fast for single residues). However, it can work (fast) only with 1-integer ids of residues.
    (or use numpy and structs, but author_insertion_code has a length limit??)
    """
    helices_start: List[int]
    helices_end: List[int]
    strands_start: List[int]
    strands_end: List[int]

    def _is_ss_for_residue(self, residue_serial: int, ss_start: List[int], ss_end: List[int]) -> bool:
        i = bisect_left(ss_start, residue_serial)
        return i != len(ss_start) and ss_end[i] >= residue_serial

    def ss_for_residue(self, residue: Residue):
        residue_serial = residue.id[1]  #  davam author_seq_id nekam, kde mam mit label_seq (nebo to predelej na author, ale vykaslat se na insertion code)

        for ss_type, ss_start, ss_end in ((SSType.HELIX, self.helices_start, self.helices_end), (SSType.STRAND, self.strands_start, self.strands_end)):
            if self._is_ss_for_residue(residue_serial, ss_start, ss_end):
                return ss_type

        return SSType.NO_TYPE


@dataclass
class SSForStructure:
    ss_for_chains: Dict[str, SSForChain]

    def ss_for_residue(self, residue: Residue):
        return self.ss_for_chains[residue.get_parent().id].ss_for_residue(residue)


class GetSecondaryStructureForStructure(CachedAnalyzer):
    """ caches SS for the whole structure

    (pdbe API response already contains the whole structure and it's not a lot of data
    however I could use endpoint also by entity id, not structure, but that has disadvantages -- To begin with, using BioPython's mmcif parser
    I don't have the entity id.)
    """

    def run(self, pdb_code: str) -> SSForStructure:
        ss_for_chains: Dict[str, SSForChain] = {}

        for molecule in get_secondary_structure(pdb_code):
            for chain in molecule['chains']:
                helices_start = []
                helices_end = []
                for helix_segment in chain['secondary_structure']['helices']:
                    helices_start.append(helix_segment['start']['residue_number'])  # residue number is WRONG, I only have author_residue_number, so either get res. number from mmcif or use auth.
                    helices_end.append(helix_segment['end']['residue_number'])

                strands_start = []
                strands_end = []
                for helix_segment in chain['secondary_structure']['strands']:
                    strands_start.append(helix_segment['start']['residue_number'])  # residue number is WRONG, I only have author_residue_number, so either get res. number from mmcif or use auth.
                    strands_end.append(helix_segment['end']['residue_number'])
                    # todo sheet_id important?

                ss_for_chains[chain['chain_id']] = SSForChain(sorted(helices_start), sorted(helices_end), sorted(strands_start), sorted(strands_end))
                # don't know if the segment order from api guaranteed ascending, so sorted

        return SSForStructure(ss_for_chains)


class GetDomainsForStructure(CachedAnalyzer):
    """ caches Domain mappings for the whole structure """

    def run(self, pdb_code: str) -> SSForStructure:
        ss_for_chains: Dict[str, SSForChain] = {}

        for molecule in get_secondary_structure(pdb_code):
            for chain in molecule['chains']:
                helices_start = []
                helices_end = []
                for helix_segment in chain['secondary_structure']['helices']:
                    helices_start.append(helix_segment['start']['residue_number'])  # residue number is WRONG, I only have author_residue_number, so either get res. number from mmcif or use auth.
                    helices_end.append(helix_segment['end']['residue_number'])

                strands_start = []
                strands_end = []
                for helix_segment in chain['secondary_structure']['strands']:
                    strands_start.append(helix_segment['start']['residue_number'])  # residue number is WRONG, I only have author_residue_number, so either get res. number from mmcif or use auth.
                    strands_end.append(helix_segment['end']['residue_number'])
                    # todo sheet_id important?

                ss_for_chains[chain['chain_id']] = SSForChain(sorted(helices_start), sorted(helices_end), sorted(strands_start), sorted(strands_end))
                # don't know if the segment order from api guaranteed ascending, so sorted

        return SSForStructure(ss_for_chains)


ResidueId = Tuple[int, str]


# TResidueData = TypeVar('TResidueData')
@dataclass(frozen=True, eq=False)
class SetOfResidues:
    # TResidueId = Tuple[int, str]  # author_seq_id, author_insertion_code? As in BioPython, but no hetero flag

    data: List[Residue]
    structure_id: str

    def get_full_id(self):
        raise NotImplementedError

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def serialize(self):
        return self.get_full_id()

    def get_atoms(self):
        """ for compatibility with FreeSASA BioPython binding, which uses get_atoms() """
        for residue in self.data:
            yield from residue

    def __add__(self, other: 'SetOfResidues'):
        return CombinedSetOfResidues.create_from_groups_of_residues(self, other)

    def __hash__(self):
        return hash(self.get_full_id())  # zatim asi hack, do budoucna trochu jinak bude cachovani fungovat (tohle by stejně asi nefunguvalo -- kdyby tam sla instance, nemohl by to GC odstranit

    def __eq__(self, other):
        return self.get_full_id() == other.get_full_id()


@dataclass(frozen=True, eq=False)
class CombinedSetOfResidues(SetOfResidues):
    combined_id: Tuple[Any]

    @classmethod
    def create_from_groups_of_residues(cls, *groups_of_residues: SetOfResidues):
        all_residues = []

        for residues in groups_of_residues:
            all_residues.extend(residues)

        return CombinedSetOfResidues(all_residues,
                                     groups_of_residues[0].structure_id,
                                     tuple(g.get_full_id() for g in groups_of_residues))

    def get_full_id(self):
        return self.combined_id


@dataclass(frozen=True, eq=False)
class ChainResidues(SetOfResidues):
    chain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id)


@dataclass(frozen=True, eq=False)
class DomainResidues(SetOfResidues):
    chain_id: str
    domain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id, self.domain_id)


class GetSASAForStructure(CachedAnalyzer):
    """ Return solvent accessible surface area for a group of residues.

    Uses FreeSASA. TODO - outputs a lot of warnings: FreeSASA: warning: guessing that atom 'CG' is symbol ' C'
    - either provide my own classifier, or temporarily redirect stderr (but then I lose error reporting)
    """

    def run(self, residues: SetOfResidues) -> float:
        sasa_structure = freesasa.structureFromBioPDB(residues)
        result = freesasa.calc(sasa_structure)
        return result.totalArea()


class GetInterdomainSurface(Analyzer):
    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_sasa: GetSASAForStructure) -> float:
        return 1/2 * (get_sasa(residues1) + get_sasa(residues2) - get_sasa(residues1 + residues2))


class CompareSecondaryStructure(SerializableCachedAnalyzer):
    """ caches the API response for the whole structure, as returned by `get_secondary_structure` """

    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_ss: GetSecondaryStructureForStructure) -> float:
        """ sequences of residues have to correspond to each other

         :returns ratio of residues with the same SS, if no SS is known for a residue, it is assumed to have NO_TYPE SS.
        """

        assert len(residues1) == len(residues2)

        ss1 = get_ss(residues1.structure_id)
        ss2 = get_ss(residues2.structure_id)

        same_count = 0
        for r1, r2 in zip(residues1, residues2):
            if ss1.ss_for_residue(r1) == ss2.ss_for_residue(r2):
                same_count += 1

        return same_count / len(residues1)


def get_c_alpha_coords(residues):
    return np.array([res['CA'].get_coord() for res in residues])


class GetRMSD(SerializableCachedAnalyzer):
    """nebo GetRMSD - naming. todo bude brát nejspíš numpy array, nebo nějakej jeho wrapper"""

    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_main_chain) -> float:
        P, Q = map(get_c_alpha_coords, (residues1, residues2))

        return rmsd.kabsch_rmsd(P, Q)
