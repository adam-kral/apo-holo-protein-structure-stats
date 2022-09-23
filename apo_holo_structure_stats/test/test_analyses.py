import logging
from pathlib import Path
from unittest import TestCase

from Bio.PDB import MMCIFParser, is_aa
import numpy as np
from scipy.spatial.transform.rotation import Rotation

from apo_holo_structure_stats.core.analyses import GetSASAForStructure, GetInterfaceBuriedArea, GetRMSD, \
    GetCenteredCAlphaCoords, GetCAlphaCoords, GetCentroid, GetRotationMatrix, GetHingeAngle
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds
from apo_holo_structure_stats.core.dataclasses import ChainResidues, DomainResidueMapping, DomainResidues, SetOfResidues

# debug aligner for preliminary sequence analysis (apo-holo/holo-holo), to see how they differ, if they differ
from Bio import Align

from input.download import parse_mmcif

aligner = Align.PairwiseAligner(mode='global',
                                open_gap_score=-0.5,
                                extend_gap_score=-0.1,
                                end_gap_score=0,)  # match by default 1 and mismatch 0


def sequences_same(rs1: SetOfResidues, rs2: SetOfResidues) -> bool:
    seq1 = [r.resname for r in rs1]
    seq2 = [r.resname for r in rs2]

    if seq1 != seq2:
        # debug print to see how seqs differ

        alignment = next(aligner.align(seq1, seq2))
        logging.info('Sequences differ, alignment:')
        logging.info(alignment)
        return False

    return True


class TestAnalyses(TestCase):
    TEST_STRUCTURE_DIR = Path(__file__).parent / 'test_data'

    def setUp(self) -> None:
        parser = MMCIFParser()
        self.structure_parse_result = self.load_test_structure('2gdu')
        self.structure = self.structure_parse_result.structure
        self.counter = 0

    def get_test_structure(self):
        s = self.structure.copy()
        s.id = f'{s.id}_{self.counter}'  # assign a unique id for the structure
        self.counter += 1
        return s

    def load_test_structure(self, pdb_code):
        return parse_mmcif(pdb_code, self.TEST_STRUCTURE_DIR / f'{pdb_code}.cif', allow_download=False)

    def test_interdomain_surface_paper(self):
        s1_parse_result = self.load_test_structure('1vr6')
        s2_parse_result = self.load_test_structure('1rzm')

        s1 = s1_parse_result.structure
        s2 = s2_parse_result.structure
        s1_mapping = s1_parse_result.bio_to_mmcif_mappings
        s2_mapping = s2_parse_result.bio_to_mmcif_mappings

        # analyze just A chain (like in the paper)
        s1_chain_a = s1[0]['A']
        s2_chain_a = s2[0]['A']

        c1_mapping = s1_mapping[0]['A']
        c2_mapping = s2_mapping[0]['A']

        # divide into domains exactly like in the paper
        # todo check that label_seq_ids do correspond to those (changed it now) NO they dont, use bio residue mapping

        # todo s1d1, ....
        s1d1 = DomainResidueMapping('D1', 'A', [13], [13 - 1 + 64])  # can reuse the domain for both structures (here for testing purposes)
        s2d1 = DomainResidueMapping('D1', 'A', [1], [64])  # can reuse the domain for both structures (here for testing purposes)

        s1d2 = DomainResidueMapping('D2', 'A', [13 - 1 + 65], [13 - 1 + 338])
        s2d2 = DomainResidueMapping('D2', 'A', [65], [338])

        s1d1 = DomainResidues.from_domain(s1d1, s1[0], c1_mapping)
        s1d2 = DomainResidues.from_domain(s1d2, s1[0], c1_mapping)

        s2d1 = DomainResidues.from_domain(s2d1, s2[0], c2_mapping)
        s2d2 = DomainResidues.from_domain(s2d2, s2[0], c2_mapping)

        interdomain_surface_computer = GetInterfaceBuriedArea((GetSASAForStructure(),))
        apo__domain_interface_area = interdomain_surface_computer(s1d1, s1d2)
        holo__domain_interface_area = interdomain_surface_computer(s2d1, s2d2)

        self.assertAlmostEqual(288, apo__domain_interface_area, delta=0.3 * 288)  # 218
        self.assertAlmostEqual(1024, holo__domain_interface_area, delta=0.3 * 1024)  # 933

        # rmsd
        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

        print(get_rmsd(s1d1+s1d2, s2d1+s2d2))  # 10.1 vs 8.0 in paper, todo celkem velky rozdil...

        # hinge
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
        screw_motion = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)
        self.assertAlmostEqual(147, 180/np.pi * screw_motion.angle, delta=2)

    def test_get_sasa_for_structure(self):
        s = self.get_test_structure()
        chain_b = s[0]['B']
        residues = ChainResidues(list(chain_b), s.id, 'B')

        sasa_computer = GetSASAForStructure()
        chain_b_sasa = sasa_computer(residues)

        self.assertGreater(chain_b_sasa, 1)

    def test_get_sasa_for_structure_with_hydrogens(self):
        parse_result = self.load_test_structure('7amj')

        residues = ChainResidues.from_chain(parse_result.structure[0]['A'], parse_result.bio_to_mmcif_mappings[0]['A'])

        sasa_computer = GetSASAForStructure()
        chain_b_sasa = sasa_computer(residues)

        self.assertGreater(chain_b_sasa, 1)

    def test_rmsd_guanylate_kinase_paper(self):
        apo_parse_result = self.load_test_structure('1ex6')
        holo_parse_result = self.load_test_structure('1ex7')

        apo = apo_parse_result.structure
        holo = holo_parse_result.structure

        apo_chain = apo[0]['B']  # note that different chain (as by dyndom), why?
        holo_chain = holo[0]['A']
        apo_mapping = apo_parse_result.bio_to_mmcif_mappings[0]['B']
        holo_mapping = holo_parse_result.bio_to_mmcif_mappings[0]['A']

        # analyze just A chain (like in the paper)
        apo_residues = ChainResidues.from_chain(apo_chain, apo_mapping)
        holo_residues = ChainResidues.from_chain(holo_chain, holo_mapping)

        logging.root.setLevel(logging.INFO)
        self.assertTrue(sequences_same(apo_residues, holo_residues))

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, GetRotationMatrix((get_centered_c_alpha_coords,))))

        rmsd = get_rmsd(apo_residues, holo_residues)
        self.assertAlmostEqual(4.4, rmsd, delta=0.1)  # 4.37 (vs 4.4 Å in the paper)

    def test_hinge_guanylate_kinase_paper(self):
        apo_parse_result = self.load_test_structure('1ex6')
        holo_parse_result = self.load_test_structure('1ex7')
        apo = apo_parse_result.structure[0]
        holo = holo_parse_result.structure[0]

        apo_chain = apo['B']  # note that different chain (as by dyndom), why?
        holo_chain = holo['A']
        apo_mapping = apo_parse_result.bio_to_mmcif_mappings[0]['B']
        holo_mapping = holo_parse_result.bio_to_mmcif_mappings[0]['A']

        apo_d1 = DomainResidues.from_domain(DomainResidueMapping('D1', 'B', [1, 84], [32, 186]), apo, apo_mapping)
        apo_d2 = DomainResidues.from_domain(DomainResidueMapping('D2', 'B', [33], [83]), apo, apo_mapping)

        # zmena chainu preci nepomahala, tak kde je zakopany pes?

        holo_d1 = DomainResidues.from_domain(DomainResidueMapping('D1', 'A', [1, 84], [32, 186]), holo, holo_mapping)
        holo_d2 = DomainResidues.from_domain(DomainResidueMapping('D2', 'A', [33], [83]), holo, holo_mapping)

        logging.root.setLevel(logging.INFO)
        self.assertTrue(sequences_same(apo_d1 + apo_d2, holo_d1 + holo_d2))

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))

        screw_motion = get_hinge_angle(apo_d1, apo_d2, holo_d1, holo_d2)
        self.assertAlmostEqual(47, 180/np.pi*screw_motion.angle, delta=0.2)  # in paper: dyndom: 47°, their principal axes 43.9

    def test_rmsd_pheromone_binding_protein_paper(self):
        apo_parse_result = self.load_test_structure('2fjy')
        holo_parse_result = self.load_test_structure('1dqe')

        # analyze just A chain (like in the paper)
        apo = apo_parse_result.structure
        holo = holo_parse_result.structure
        apo_chain = apo[0]['A']
        holo_chain = holo[0]['A']
        apo_mapping = apo_parse_result.bio_to_mmcif_mappings[0]['A']
        holo_mapping = holo_parse_result.bio_to_mmcif_mappings[0]['A']

        # # logging.root.setLevel(logging.INFO)
        # self.assertTrue(sequences_same(apo_chain, holo_chain))

        apo_residues = ChainResidues.from_chain(apo_chain, apo_mapping)
        holo_residues = ChainResidues.from_chain(holo_chain, holo_mapping)

        # chains also have leading and trailing residues not present in the other, remove them
        apo_residues = ChainResidues(apo_residues.data[:-5], apo.id, apo_chain.id)
        holo_residues = ChainResidues(holo_residues.data[6:], holo.id, holo_chain.id)

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, GetRotationMatrix((get_centered_c_alpha_coords,))))

        rmsd = get_rmsd(apo_residues, holo_residues)
        self.assertAlmostEqual(7, rmsd, delta=1)  # todo 6.24 (vs 7.0 in the paper)

    def test_interdomain_surface(self):
        s = self.get_test_structure()
        chain_a = ChainResidues(list(s[0]['A']), s.id, 'A')
        chain_b = ChainResidues(list(s[0]['B']), s.id, 'B')

        interdomain_surface_computer = GetInterfaceBuriedArea((GetSASAForStructure(),))
        area = interdomain_surface_computer(chain_a, chain_b)
        self.assertGreater(area, 1)

    def test_rmsd_translated(self):
        s = self.get_test_structure()
        chain_a_copy = s[0]['A'].copy()

        # move the copy by 1 angstrom
        for atom in chain_a_copy.get_atoms():
            atom.coord += (1, 0, 0)

        chain_a = ChainResidues([r for r in s[0]['A'] if is_aa(r)], s.id, 'A')
        chain_a_copy = ChainResidues([r for r in chain_a_copy if is_aa(r)], f'moved_{s.id}', 'A')

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, GetRotationMatrix((get_centered_c_alpha_coords,))))

        rmsd = get_rmsd(chain_a, chain_a_copy)
        self.assertAlmostEqual(0, rmsd, places=5)

    def test_rmsd_rotated_and_translated(self):
        s1 = self.get_test_structure()
        s2 = self.get_test_structure()

        # rotate and translate s2
        AXIS_DIRECTION = np.array([11,2,0])
        AXIS_DIRECTION = AXIS_DIRECTION / np.linalg.norm(AXIS_DIRECTION)  # following code expects a unit vector
        ANGLE = np.pi/4
        TRANSLATION = np.array([1, 200, 7])

        atoms = list(s2.get_atoms())
        coords = np.array([a.get_coord() for a in atoms])

        rotation = Rotation.from_rotvec(ANGLE * AXIS_DIRECTION)

        for atom, new_coord in zip(atoms, rotation.apply(coords) + TRANSLATION):
            atom.set_coord(new_coord)

        chain_a_mapping = self.structure_parse_result.bio_to_mmcif_mappings[0]['A']

        # test rmsd is still 0
        chain_a = ChainResidues.from_chain(s1[0]['A'], chain_a_mapping)
        chain_a_rotated = ChainResidues.from_chain(s2[0]['A'], chain_a_mapping)

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, GetRotationMatrix((get_centered_c_alpha_coords,))))

        rmsd = get_rmsd(chain_a, chain_a_rotated)
        self.assertAlmostEqual(0, rmsd, places=4)

    def test_get_hinge_angle(self):
        s1 = self.get_test_structure()
        s2 = self.get_test_structure()

        s2.id = f'{s1.id}_with_rotated_chain'

        chain_a_mapping = self.structure_parse_result.bio_to_mmcif_mappings[0]['A']
        chain_b_mapping = self.structure_parse_result.bio_to_mmcif_mappings[0]['B']
        # usually the domains are from the same chain, but this is just a test..
        s1d1 = ChainResidues.from_chain(s1[0]['A'], chain_a_mapping)
        s1d2 = ChainResidues.from_chain(s1[0]['B'], chain_b_mapping)
        s2d1 = ChainResidues.from_chain(s2[0]['A'], chain_a_mapping)
        s2d2 = ChainResidues.from_chain(s2[0]['B'], chain_b_mapping)

        # rotate second domain over a defined screw axis, then check if GetHingeAngle indeed computes the correct parameters (angle, translation)

        # define the screw axis
        AXIS_DIRECTION = np.array([1,2,0])
        AXIS_DIRECTION = AXIS_DIRECTION / np.linalg.norm(AXIS_DIRECTION)  # following code expects a unit vector
        AXIS_LOCATION = np.array([52.71183395385742, 44.92530822753906, -11.425999641418457])  # a random pivot point (the axis goes through it)
        ANGLE = np.pi/4
        TRANSLATION_IN_AXIS = 3

        # move along the screw axis using scipy
        s2d2_atoms = [atom for residue in s2d2 for atom in residue]
        s2d2_atom_coords = np.array([atom.coord for atom in s2d2_atoms])

        rotation = Rotation.from_rotvec(ANGLE * AXIS_DIRECTION)
        rotated_s2d2_atom_coords = rotation.apply(s2d2_atom_coords - AXIS_LOCATION) + AXIS_LOCATION

        for atom, new_coord in zip(s2d2_atoms, rotated_s2d2_atom_coords + AXIS_DIRECTION*TRANSLATION_IN_AXIS):
            atom.set_coord(new_coord)

        # compute the hinge angle with GetHingeAngle
        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, GetRotationMatrix((get_centered_c_alpha_coords,))))

        screw_motion = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)

        self.assertAlmostEqual(ANGLE, screw_motion.angle, places=3)
        self.assertAlmostEqual(TRANSLATION_IN_AXIS, screw_motion.translation_in_axis, places=3)
        # GetHingeAngle does not return AXIS_DIRECTION, or AXIS_LOCATION (yet)

    def test_get_hinge_angle_different_orientation(self):
        """ Orient the two structures differently. """
        s1 = self.get_test_structure()

        # orient the second structure arbitrarily, but differently that s1 (for testing purposes)
        s2 = self.get_test_structure()
        s2.id = f'{s1.id}_with_rotated_chain'

        ANGLE = np.pi * 0.666
        AXIS_DIRECTION = np.array([4, 7, 13])
        AXIS_DIRECTION = AXIS_DIRECTION / np.linalg.norm(AXIS_DIRECTION)  # following code expects a unit vector
        rotation = Rotation.from_rotvec(ANGLE * AXIS_DIRECTION)
        for atom in s2.get_atoms():
            atom.set_coord(rotation.apply(atom.get_coord()) + [6, 6, 6])

        chain_a_mapping = self.structure_parse_result.bio_to_mmcif_mappings[0]['A']
        chain_b_mapping = self.structure_parse_result.bio_to_mmcif_mappings[0]['B']
        # usually the domains are from the same chain, but this is just a test..
        s1d1 = ChainResidues.from_chain(s1[0]['A'], chain_a_mapping)
        s1d2 = ChainResidues.from_chain(s1[0]['B'], chain_b_mapping)
        s2d1 = ChainResidues.from_chain(s2[0]['A'], chain_a_mapping)
        s2d2 = ChainResidues.from_chain(s2[0]['B'], chain_b_mapping)

        # rotate second domain over a defined screw axis, then check if GetHingeAngle indeed computes the correct parameters (angle, translation)

        # define the screw axis
        AXIS_DIRECTION = np.array([1,2,0])
        AXIS_DIRECTION = AXIS_DIRECTION / np.linalg.norm(AXIS_DIRECTION)  # following code expects a unit vector
        AXIS_LOCATION = np.array([52.71183395385742, 44.92530822753906, -11.425999641418457])  # a random pivot point (the axis goes through it)
        ANGLE = np.pi/4
        TRANSLATION_IN_AXIS = 3

        # move along the screw axis using scipy
        s2d2_atoms = [atom for residue in s2d2 for atom in residue]
        s2d2_atom_coords = np.array([atom.coord for atom in s2d2_atoms])

        rotation = Rotation.from_rotvec(ANGLE * AXIS_DIRECTION)
        rotated_s2d2_atom_coords = rotation.apply(s2d2_atom_coords - AXIS_LOCATION) + AXIS_LOCATION

        for atom, new_coord in zip(s2d2_atoms, rotated_s2d2_atom_coords + AXIS_DIRECTION*TRANSLATION_IN_AXIS):
            atom.set_coord(new_coord)

        # compute the hinge angle with GetHingeAngle
        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, GetRotationMatrix((get_centered_c_alpha_coords,))))

        screw_motion = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)

        self.assertAlmostEqual(ANGLE, screw_motion.angle, places=3)
        self.assertAlmostEqual(TRANSLATION_IN_AXIS, screw_motion.translation_in_axis, places=3)
        # GetHingeAngle does not return AXIS_DIRECTION, or AXIS_LOCATION (yet)
