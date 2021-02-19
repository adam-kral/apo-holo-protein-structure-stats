from bisect import bisect_left, bisect
from enum import Enum, auto
import itertools
from dataclasses import dataclass
from typing import Iterator, List, Dict, Tuple, Any, DefaultDict, NamedTuple, TypeVar, Generic

import freesasa as freesasa
import numpy as np
import rmsd
from Bio.PDB import is_aa, NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

from apo_holo_structure_stats.input.download import get_secondary_structure, get_domains
from .base_analyses import CachedAnalyzer, Analyzer, SerializableCachedAnalyzer, SerializableAnalyzer


freesasa.setVerbosity(freesasa.nowarnings)  # FreeSASA: warning: guessing that atom 'CB' is symbol ' C' ..., or todo can set a custom classifier?


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


class ResidueId(NamedTuple):
    """ Allows unique identification of a residue within a PDB structure

     Esp. within Bio.PDB.Structure (otherwise it'd be better with label_seq_id and label_asym_id) """
    auth_seq_id: int
    insertion_code: str

    chain_id: str  # is already present in e.g. ChainResidueData, but usually I expect just SetOfResidueData (may be from multiple chains,
    # simplest, ok)

    @classmethod
    def from_bio_residue(cls, residue: Residue, chain_id: str):
        return ResidueId(residue.id[1], residue.id[2], chain_id)


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

    def ss_for_residue(self, r: ResidueId):
        residue_serial = r.auth_seq_id  #  davam author_seq_id nekam, kde mam mit label_seq (nebo to predelej na author, ale vykaslat se na insertion code)

        for ss_type, ss_start, ss_end in ((SSType.HELIX, self.helices_start, self.helices_end), (SSType.STRAND, self.strands_start, self.strands_end)):
            if self._is_ss_for_residue(residue_serial, ss_start, ss_end):
                return ss_type

        return SSType.NO_TYPE


@dataclass
class SSForStructure:
    ss_for_chains: Dict[str, SSForChain]

    def ss_for_residue(self, r: ResidueId):
        return self.ss_for_chains[r.chain_id].ss_for_residue(r)


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
                    helices_start.append(helix_segment['start']['author_residue_number'])  # or 'residue_number', but I don't have label_seq_id in BioPython-parsed Structures, Also I ignore 'author_insertion_code'
                    helices_end.append(helix_segment['end']['author_residue_number'])

                strands_start = []
                strands_end = []
                for helix_segment in chain['secondary_structure']['strands']:
                    strands_start.append(helix_segment['start']['author_residue_number'])
                    strands_end.append(helix_segment['end']['author_residue_number'])
                    # todo sheet_id important?

                ss_for_chains[chain['chain_id']] = SSForChain(sorted(helices_start), sorted(helices_end), sorted(strands_start), sorted(strands_end))
                # don't know if the segment order from api guaranteed ascending, so sorted

        return SSForStructure(ss_for_chains)




TResidueData = TypeVar('TResidueData')

@dataclass(frozen=True, eq=False)
class SetOfResidueData(Generic[TResidueData]):
    # TResidueId = Tuple[int, str]  # author_seq_id, author_insertion_code? As in BioPython, but no hetero flag

    data: List[TResidueData]
    structure_id: str

    def get_full_id(self):
        raise NotImplementedError

    def __len__(self):
        return len(self.data)

    def __iter__(self) -> Iterator[TResidueData]:
        return iter(self.data)

    def serialize(self):
        return self.get_full_id()

    def __add__(self, other: 'SetOfResidueData[TResidueData]'):
        return CombinedSetOfResidueData[TResidueData].create_from_groups_of_residues(self, other)

    def __hash__(self):
        return hash(self.get_full_id())  # zatim asi hack, do budoucna trochu jinak bude cachovani fungovat (tohle by stejně asi nefunguvalo -- kdyby tam sla instance, nemohl by to GC odstranit

    def __eq__(self, other):
        return self.get_full_id() == other.get_full_id()


@dataclass(frozen=True, eq=False)
class CombinedSetOfResidueData(SetOfResidueData[TResidueData]):
    combined_id: Tuple[Any]

    @classmethod
    def create_from_groups_of_residues(cls, *groups_of_residues: SetOfResidueData[TResidueData]):
        all_residues = []

        for residues in groups_of_residues:
            all_residues.extend(residues)

        return CombinedSetOfResidueData[TResidueData](all_residues,
                                     groups_of_residues[0].structure_id,
                                     tuple(sorted(g.get_full_id() for g in groups_of_residues)))  # sorted, so __eq__ works as expected

    def get_full_id(self):
        return self.combined_id


@dataclass(frozen=True, eq=False)
class ChainResidueData(SetOfResidueData[TResidueData]):
    chain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id)


@dataclass(frozen=True, eq=False)
class DomainResidueData(SetOfResidueData[TResidueData]):
    chain_id: str
    domain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id, self.domain_id)


SetOfResidues = SetOfResidueData[Residue]


class ChainResidues(ChainResidueData[Residue]):
    @classmethod
    def from_bio_chain(cls, chain: Chain):
        return cls([r for r in chain if is_aa(r)], chain.get_parent().get_parent().id, chain.id)


@dataclass
class DomainResidueMapping:
    domain_id: str
    chain_id: str
    segment_beginnings: List[int]
    segment_ends: List[int]

    def __iter__(self) -> Iterator[int]:
        """ Returns a sequence of residue auth_seq_id in a domain """
        for segment_start, segment_end in zip(self.segment_beginnings, self.segment_ends):
            yield from range(segment_start, segment_end + 1)

    def __contains__(self, auth_seq_id):
        domain_segment_index = -1 + bisect(self.segment_beginnings, auth_seq_id)
        return domain_segment_index >= 0 and auth_seq_id <= self.segment_ends[domain_segment_index]

    def __len__(self):
        residue_count = 0
        for segment_start, segment_end in zip(self.segment_beginnings, self.segment_ends):
            residue_count += segment_end + 1 - segment_start

        return residue_count

    def to_set_of_residue_ids(self, structure_id: str) -> SetOfResidueData[ResidueId]:
        return DomainResidueData[ResidueId]([ResidueId(auth_seq_id, ' ', self.chain_id) for auth_seq_id in self], structure_id, self.chain_id, self.domain_id)


class GetDomainsForStructure(CachedAnalyzer):
    """ caches Domain mappings for the whole structure """

    def run(self, pdb_code: str) -> List[DomainResidueMapping]:
        domains = DefaultDict[str, DomainResidueMapping]()

        for superfamily_id, superfamily in get_domains(pdb_code).items():
            for domain_segment in superfamily['mappings']:
                domain_id = domain_segment['domain']

                if domain_id not in domains:
                    domains[domain_id] = DomainResidueMapping(domain_id, domain_segment['chain_id'], [], [])

                domains[domain_id].segment_beginnings.append(domain_segment['start']['author_residue_number'])
                domains[domain_id].segment_ends.append(domain_segment['end']['author_residue_number'])

        return list(domains.values())


class DomainResidues(DomainResidueData[Residue]):
    @classmethod
    def from_domain(cls, domain: DomainResidueMapping, bio_structure: Structure):
        bio_chain = bio_structure[0][domain.chain_id]

        try:
            domain_residues = [bio_chain[' ', auth_seq_id, ' '] for auth_seq_id in domain]
        except KeyError:
            # a residue probably has a hetero flag (!= ' '), e.g. 'H_CSU' (CYSTEINE-S-SULFONIC ACID)
            # formally this could be the only code (try not needed), I don't know which is faster
            domain_residues = [residue for residue in bio_chain if residue.get_id()[1] in domain]

        return cls(domain_residues, bio_structure.id, domain.chain_id, domain.domain_id)


class GetSASAForStructure(CachedAnalyzer):
    """ Return solvent accessible surface area for a group of residues.

    Uses FreeSASA. TODO - outputs a lot of warnings: FreeSASA: warning: guessing that atom 'CG' is symbol ' C'
    - either provide my own classifier, or disable all freesasa warnings (done at the file head)
    """

    def run(self, residues: SetOfResidues) -> float:
        # attach method get_atoms used by freesasa's BioPython binding (so that it behaves like BioPython's Entity)
        def get_atoms(self):
            for r in self:
                yield from r.get_atoms()

        bound_method =  get_atoms.__get__(residues)
        object.__setattr__(residues, 'get_atoms', bound_method)  # setting to a _frozen_ dataclass (SetOfResidues)

        # use freesasa to compute SASA
        sasa_structure = freesasa.structureFromBioPDB(residues)
        result = freesasa.calc(sasa_structure)
        return result.totalArea()


class GetInterdomainSurface(SerializableAnalyzer):
    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_sasa: GetSASAForStructure) -> float:
        return 1/2 * (get_sasa(residues1) + get_sasa(residues2) - get_sasa(residues1 + residues2))  # the values in the paper seem not to be multiplied by 1/2


class CompareSecondaryStructure(SerializableAnalyzer):
    """ caches the API response for the whole structure, as returned by `get_secondary_structure` """

    def run(self, residues1: SetOfResidueData[ResidueId], residues2: SetOfResidueData[ResidueId], get_ss: GetSecondaryStructureForStructure) -> float:
        """ sequences of residues have to correspond to each other

         :returns ratio of residues with the same SS, if no SS is known for a residue, it is assumed to have NO_TYPE SS.
        """

        assert len(residues1) == len(residues2)

        ss1 = get_ss(residues1.structure_id)
        ss2 = get_ss(residues2.structure_id)

        same_count = 0
        for r1_id, r2_id in zip(residues1, residues2):
            if ss1.ss_for_residue(r1_id) == ss2.ss_for_residue(r2_id):
                same_count += 1

        return same_count / len(residues1)


class GetCAlphaCoords(CachedAnalyzer):
    def run(self, residues: SetOfResidues) -> np.ndarray:
        return np.array([res['CA'].get_coord() for res in residues])

""" nevyhody tohoto postupu --> neda se moc composovat (volat) s parametry zjistenymi ruznymi fcemi. Napr. do GetCentroid nemuzu jednou
poslat souradnice CA podruhy napr vsech atomu, protoze GetCentroid si musi volat tu funkci get coords samo. Ale vlastne bych tam mohl poslat jinou
dependency (takze i typ, to se ale da vyresit), (ale co má stejný signature, to bych normálně nemusel)"""


class GetCentroid(CachedAnalyzer):
    def run(self, residues: SetOfResidues, get_c_alpha_coords: GetCAlphaCoords) -> float:
        return rmsd.centroid(get_c_alpha_coords(residues))


class GetCenteredCAlphaCoords(CachedAnalyzer):
    def run(self, residues: SetOfResidues, get_c_alpha_coords: GetCAlphaCoords, get_centroid: GetCentroid) -> np.ndarray:
        return get_c_alpha_coords(residues) - get_centroid(residues)


class GetRotationMatrix(CachedAnalyzer):
    """ Returns a matrix corresponding to rotation of P onto Q minimizing RMSD (the matrix is actually transposed, see usage below).

    Q fixed, P would move onto Q, if done the following: P = np.dot(P, get_rotation_matrix(P, Q))

    Rotation matrix is transposed, as the point coordinates are too (as row vectors); this is the practice of
    the rmsd package.
    """
    def run(self, r1: SetOfResidues, r2: SetOfResidues, get_centered_c_alpha_coords: GetCenteredCAlphaCoords):
        P, Q = map(get_centered_c_alpha_coords, (r1, r2))

        return rmsd.kabsch(P, Q)


class GetRMSD(SerializableCachedAnalyzer):
    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_centered_c_alpha_coords: GetCenteredCAlphaCoords, get_rotation_matrix: GetRotationMatrix) -> float:
        P, Q = map(get_centered_c_alpha_coords, (residues1, residues2))

        # following code from rmsd.kabsch_rmsd()
        P = np.dot(P, get_rotation_matrix(residues1, residues2))
        return rmsd.rmsd(P, Q)


class GetHingeAngle(SerializableCachedAnalyzer):
    """ Computes the angle of a two-domain movement. Provide two domains from one structure and the corresponding ones from the second structure.

    Uses first domain as a reference and computes an angle by which the second domain tilted to get from 1st structure -> 2nd structure.
    (Returns absolute values, so the result is the same if structures/domains are interchanged.)
    Can compute screw axis parameters for the movement. (However, currently computation of the location of the screw axis is missing.)
    """
    def run(self, s1d1: SetOfResidues, s1d2: SetOfResidues, s2d1: SetOfResidues, s2d2: SetOfResidues,
            get_c_alpha_coords: GetCAlphaCoords,
            get_centroid: GetCentroid,
            get_rotation_matrix: GetRotationMatrix) -> (float, float):
        """ :param s1d1: a domain in the first structure
        :param s1d2: another domain in the first structure
        :param s2d1: a domain in the second structure corresponding to s1d1 domain (same length required)
        :param s2d2: a domain in the second structure corresponding to s1d2 domain (same length required)
        :return: tuple of angle (non-negative -- absolute value) and translation (again non-negative).
        Angle (movement of a domain between the structures) and translation -> in the direction of the rotation axis
        """

        # superimpose s2d1 on s1d1 (corr. domains, 1st domain chosen arbitrarily), taking s2d2 along with it
        # -> translate and rotate

        # superimpose the moved s2d2 on s1d2, remember translation+rotation. Decompose translation+rotation into rotation along a single
        # axis somewhere in space + translation along that axis. (screw motion).

        # axis direction: 1 eigenvector of rot. matrix. Translation along that axis: subtract from translation it's projection
        # to screw axis _|_ space (P_s_|_).
        # [ In case we would need the location of the screw axis: the projection P_s_|_ is the length of
        # the chord line (how the centroid rotated along the positioned screw axis). I know the rot. angle, so the location of screw axis is
        # at two possible positions, above and below the chord line, however the angle sign should solve this
        # r = chord/(2 * sin(angle/2)), compute triangle height <r,chord,r> and multiply by that a unit vector perpendiculat to the chord,
        # and in the rotation plane. This will be the candidate (one of two, see the angle sign) for the screw axis location.
        # (A point lying on the axis). ]

        # superimpose the two structures by their first domain
        U = get_rotation_matrix(s2d1, s1d1)

        s1d2_coords = get_c_alpha_coords(s1d2)
        fixed_s2d2_coords = np.matmul(get_c_alpha_coords(s2d2) - get_centroid(s2d1), U) + get_centroid(s1d1)

        # compare d2 coords
        fixed_s2d2_coords_centroid = fixed_s2d2_coords.mean(axis=0)

        d2_rotation_m = rmsd.kabsch(fixed_s2d2_coords - fixed_s2d2_coords_centroid, s1d2_coords - get_centroid(s1d2))
        d2_translation = get_centroid(s1d2) - fixed_s2d2_coords_centroid

        # calculate return values
        angle = np.arccos((np.trace(d2_rotation_m) - 1) / 2)  # formula tr(Rot_matrix) = 1 + 2cos(angle)

        rotation_axis_plane = d2_rotation_m - np.eye(3)  # rotation axis is the normal of that plane, row/col vectors generate the plane
        axis = np.cross(rotation_axis_plane[:, 0], rotation_axis_plane[:, 1])  # a vector in the direction of the axis
        # (result of cross product is perpendicular to both vectors, these are in rotation_axis_plane) https://math.stackexchange.com/q/2074316/331963
        axis = axis / np.linalg.norm(axis)  # make it a unit vector

        # translation projected into the rotation axis (get the translation component of a screw motion)
        translation_in_axis = np.dot(d2_translation, axis)

        # it could also be possible to calculate here the location of rotation axis --> so we would have the complete screw description
        # of the domain movement, see [ ] above

        return abs(angle), abs(translation_in_axis)
