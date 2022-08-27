import itertools
from typing import Iterator, List, Dict

import freesasa as freesasa
import numpy as np
import rmsd
from Bio.PDB import is_aa, NeighborSearch
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from apo_holo_structure_stats.input.download import get_secondary_structure, get_domains
from settings import Settings
from .base_analyses import CachedAnalyzer, SerializableCachedAnalyzer, SerializableAnalyzer
from .dataclasses import SSForChain, SSForStructure, SetOfResidueData, SetOfResidues, DomainResidueMapping, \
    ScrewMotionResult, ChainResidues
from .biopython_to_mmcif import ResidueId

freesasa.setVerbosity(freesasa.nowarnings)  # FreeSASA: warning: guessing that atom 'CB' is symbol ' C' ..., or todo can set a custom classifier?


class AnalysisException(Exception):
    pass


class MissingDataException(AnalysisException):
    """ Can be raised for example as __cause__`d by APIException. """
    pass


def get_hetero_atom_residues(struct: Entity) -> Iterator[Residue]:
    """ non-polymer ligands, excl. water

    BioPython's residue is a group of atoms belonging to the same molecule (made possible by mmcifs author_seq_id of heteroatoms,
    label_seq_id is '.' because a heterocompound does not belong to the poly-chain)
    scans all chains (we deliberately ignore author chain assignment in case of a heteroatom)
    # todo maybe this is wrong, as sub-compound are covalently bound into a single ligand, if they share a chain assignment?
        seems true for saccharides (but they do have monomers), they share label_asym_id, but not the auth_asym_id... (as in biopython)
        in the paper, they mark multimeric with '-' between monomers and not covalently bound, distinct, ligands probably delimited
        with ';'
        """

    return filter(lambda res: res.id[0].startswith('H_'), struct.get_residues())  # .id[0] is a residue's hetatm flag


def get_short_peptide_ligands(struct: Entity, peptide_length_limit: int) -> Iterator[ChainResidues]:
    peptide_chains = filter(lambda chain: sum(is_aa(residue) for residue in chain) <= peptide_length_limit, struct.get_chains())
    return [ChainResidues([r for r in peptide_chain if is_aa(r)],
                          peptide_chain.get_parent().get_parent().id,
                          peptide_chain.id) for peptide_chain in peptide_chains]


def get_all_ligands(struct: Model) -> Iterator[Entity]:
    return itertools.chain(
        get_hetero_atom_residues(struct),
        get_short_peptide_ligands(struct, Settings.LigandSpec.PEPTIDE_MAX_LENGTH)
    )


def get_defined_ligands(struct: Model, chain: SetOfResidues) -> Iterator[Entity]:
    """ Get defined ligands bound with defined specificity.  (See Settings.LigandSpec)

    :param struct: The whole resolved structure, ligands are identified there.
    :param chain: Only the polypeptide chain (no HETATMs etc.), as the specificity is calculated wrt. that.
    :return: ligands (Either Bio.Residue or ChainResidues in case of a peptide ligand)
    """
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
    ligands = list(filter(lambda lig: has_at_least_n_non_hydrogen_atoms(lig, Settings.LigandSpec.MIN_NON_H_ATOMS),
                          get_all_ligands(struct)))

    # ligand is within RADIUS in contact with MIN_RESIDUES_WITHIN_LIGAND residues

    # (in original paper they used a program LPC, for ensuring specific interaction of ligand with at least 6 residue,
    # this is a "shortcut", a temporary condition (simple))

    chain_atoms = list(chain.get_atoms())
    ns = NeighborSearch(chain_atoms)

    MIN_RESIDUES_WITHIN_LIGAND = Settings.LigandSpec.MIN_RESIDUES_WITHIN_LIGAND
    RADIUS = Settings.LigandSpec.MIN_RESIDUES_WITHIN_LIGAND__RADIUS
    # todo calculate average number of protein heavy atoms in 4.5 Å within ligand _atom_ (paper says 6)

    for ligand in ligands:
        residues_in_contact_with_ligand = set()  # including the ligand itself (in biopython, non-peptide ligand is
        # in the same chain usually, but in a different residue)

        ligand_residues = set()  # residues that compose the ligand

        for ligand_atom in ligand.get_atoms():  # ligand can be a chain or a residue
            ligand_residues.add(ligand_atom.get_parent())
            chain_atoms_in_contact = ns.search(ligand_atom.get_coord(), RADIUS)

            for atom in chain_atoms_in_contact:
                # exclude hydrogen atoms (as in the paper)
                if atom.element == 'H':
                    continue

                residues_in_contact_with_ligand.add(atom.get_parent())

        # exclude the ligand itself from the set of contact residues
        residues_in_contact_with_ligand -= ligand_residues

        if len(residues_in_contact_with_ligand) >= MIN_RESIDUES_WITHIN_LIGAND:
            yield ligand


class GetChains(CachedAnalyzer):
    """ Use minimum number of observed residues (in bio.Chain) to define a chain for the subsequent analyses. """
    def run(self, struct: Model) -> List[Chain]:
        # (or could use all residues (incl. not observed) -> entity_poly_seq in BiopythonToMmcif)
        return list(filter(
            lambda chain: sum(is_aa(residue) for residue in chain) >= Settings.MIN_OBSERVED_RESIDUES_FOR_CHAIN,
            struct.get_chains()
        ))


class GetMainChain(CachedAnalyzer):
    def run(self, struct: Model, get_chains: GetChains) -> Chain:
        """ :returns one chain (the first of GetChains), because hopefully we're working with single-chain structures """
        return get_chains(struct)[0]


class IsHolo(CachedAnalyzer):
    def run(self, struct: Model, chain: SetOfResidues):
        acceptable_ligands = list(get_defined_ligands(struct, chain))
        return len(acceptable_ligands) > 0


def is_sorted(lst: List):
    return all(x <= y for x, y in zip(lst, lst[1:]))


def sort_bound_lists(*lists, by=0):
    """ Sort lists as if elements were bound across lists at each index. Sort by `by`th list.
    """
    assert 0 <= by < len(lists)
    tuples = zip(*lists)

    output_lists = [[] for _ in range(len(lists))]

    for tuple in sorted(tuples, key=lambda t: t[by]):
        for i in range(len(lists)):
            output_lists[i].append(tuple[i])

    return output_lists


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

                if 'helices' in chain['secondary_structure']:  # otherwise got error (with strands below), seems (some) fields are optional
                    for helix_segment in chain['secondary_structure']['helices']:
                        helices_start.append(helix_segment['start']['residue_number'])  # (is label_seq_id)
                        helices_end.append(helix_segment['end']['residue_number'])

                strands_start = []
                strands_end = []

                if 'strands' in chain['secondary_structure']:
                    for helix_segment in chain['secondary_structure']['strands']:
                        strands_start.append(helix_segment['start']['residue_number'])
                        strands_end.append(helix_segment['end']['residue_number'])
                        # todo sheet_id important?

                # segment order from api not guaranteed ascending, and sometimes is NOT
                helices_start, helices_end = sort_bound_lists(helices_start, helices_end)
                strands_start, strands_end = sort_bound_lists(strands_start, strands_end)

                ss_for_chains[chain['chain_id']] = SSForChain(helices_start, helices_end, strands_start, strands_end)

        return SSForStructure(ss_for_chains)


class GetDomainsForStructure(CachedAnalyzer):
    """ caches Domain mappings for the whole structure """

    def run(self, pdb_code: str) -> List[DomainResidueMapping]:
        domains: Dict[str, DomainResidueMapping] = {}

        for superfamily_id, superfamily in get_domains(pdb_code).items():
            for domain_segment in superfamily['mappings']:
                domain_id = domain_segment['domain']

                if domain_id not in domains:
                    domains[domain_id] = DomainResidueMapping(domain_id, domain_segment['chain_id'], [], [])

                domains[domain_id].segment_beginnings.append(domain_segment['start']['residue_number'])
                domains[domain_id].segment_ends.append(domain_segment['end']['residue_number'])

        # segment order from api not guaranteed ascending
        for domain in domains.values():
            domain.segment_beginnings, domain.segment_ends = sort_bound_lists(domain.segment_beginnings, domain.segment_ends)

        return list(domains.values())


def get_sasa(residues: SetOfResidues) -> freesasa.Result:
    # use freesasa to compute SASA
    # it uses method get_atoms (which is defined on SetOfResidues). It does not include hydrogen atoms, because
    # otherwise it raises an error (hydrogens should be ignored by default in freesasa anyway, but they say
    # structureFromBioPDB is experimental)
    sasa_structure = freesasa.structureFromBioPDB(residues)
    return freesasa.calc(sasa_structure)


class GetSASAForStructure(CachedAnalyzer):
    """ Return solvent accessible surface area for a group of residues.

    Uses FreeSASA. TODO - outputs a lot of warnings: FreeSASA: warning: guessing that atom 'CG' is symbol ' C'
    - either provide my own classifier, or disable all freesasa warnings (done at the file head)
    """

    def run(self, residues: SetOfResidues) -> float:
        freesasa_result = get_sasa(residues)
        return freesasa_result.totalArea()


class GetInterfaceBuriedArea(SerializableAnalyzer):
    """ Returns the area of interface between two domains (or sets of residues). In Angstroms^2. """
    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_sasa: GetSASAForStructure) -> float:
        return get_sasa(residues1) + get_sasa(residues2) - get_sasa(residues1 + residues2)


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
        ca_coords = []

        for res in residues:
            try:
                ca_coords.append(res['CA'].get_coord())
            except KeyError as e:
                raise AnalysisException(f'Residue {res.id} is missing c_alpha coordinates')  # maybe should be also
                # ValueError somehow?

        return np.array(ca_coords)

""" nevyhody tohoto postupu --> neda se moc composovat (volat) s parametry zjistenymi ruznymi fcemi. Napr. do GetCentroid nemuzu jednou
poslat souradnice CA podruhy napr vsech atomu, protoze GetCentroid si musi volat tu funkci get coords samo. Ale vlastne bych tam mohl poslat jinou
dependency (takze i typ, to se ale da vyresit), (ale co má stejný signature (interface v oop), to bych normálně nemusel)
- takže když bude čas změnit? Takhle jsem to dělal kvůli cachovaní nebo serializaci výsledků, ale cachování je zbytečný
- ukazuje se, že nejvíc času trvá načtení struktury z mmcifu a pak počítání SASA (tam bych to využil, ale marginálně). 
Zbytek je skoro nic. I protože jsem udělal to prefetch (1struct_analyses) a cachování API response get_ss a get_domains,
 to ale asi díky tomuhle právě ale.."""


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


class GetRMSD(SerializableAnalyzer):
    def run(self, residues1: SetOfResidues, residues2: SetOfResidues, get_centered_c_alpha_coords: GetCenteredCAlphaCoords, get_rotation_matrix: GetRotationMatrix) -> float:
        P, Q = map(get_centered_c_alpha_coords, (residues1, residues2))

        # following code from rmsd.kabsch_rmsd()
        P = np.dot(P, get_rotation_matrix(residues1, residues2))
        return rmsd.rmsd(P, Q)


class GetHingeAngle(SerializableAnalyzer):
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

        # superimpose s2d1 on s1d1 (correspoding domains, 1st domain chosen arbitrarily), taking s2d2 along with it
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

        return ScrewMotionResult(abs(angle), abs(translation_in_axis), float(np.linalg.norm(d2_translation)))
