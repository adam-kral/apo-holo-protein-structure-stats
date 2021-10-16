import logging
from pathlib import Path
from unittest import TestCase

from Bio.PDB import MMCIFParser, is_aa
import numpy as np
from scipy.spatial.transform.rotation import Rotation

from apo_holo_structure_stats.core.analyses import GetSASAForStructure, GetInterdomainSurface, GetRMSD, \
    GetCenteredCAlphaCoords, GetCAlphaCoords, GetCentroid, GetRotationMatrix, GetHingeAngle, DomainResidues
from apo_holo_structure_stats.core.dataclasses import ChainResidues, DomainResidueMapping
from apo_holo_structure_stats.pipeline.run_analyses import sequences_same, chain_to_polypeptide


class TestAnalyses(TestCase):
    TEST_STRUCTURE_DIR = Path(__file__).parent / 'test_data'

    def setUp(self) -> None:
        self.structure = MMCIFParser().get_structure('2gdu', self.TEST_STRUCTURE_DIR / '2gdu.cif')
        self.counter = 0

    def get_test_structure(self):
        s = self.structure.copy()
        s.id = f'{s.id}_{self.counter}'  # assign a unique id for the structure
        self.counter += 1
        return s

    def load_test_structure(self, pdb_code):
        return MMCIFParser().get_structure(pdb_code, self.TEST_STRUCTURE_DIR / f'{pdb_code}.cif')

    def test_interdomain_surface_paper(self):
        s1 = self.load_test_structure('1vr6')
        s2 = self.load_test_structure('1rzm')

        # analyze just A chain (like in the paper)
        s1_chain_a = s1[0]['A']
        s2_chain_a = s2[0]['A']

        # divide into domains exactly like in the paper
        d1 = DomainResidueMapping('D1', 'A', [1], [64])  # can reuse the domain for both structures (here for testing purposes)
        d2 = DomainResidueMapping('D2', 'A', [65], [338])

        s1d1 = DomainResidues.from_domain(d1, s1)
        s1d2 = DomainResidues.from_domain(d2, s1)

        s2d1 = DomainResidues.from_domain(d1, s2)
        s2d2 = DomainResidues.from_domain(d2, s2)

        interdomain_surface_computer = GetInterdomainSurface((GetSASAForStructure(),))
        apo__domain_interface_area = interdomain_surface_computer(s1d1, s1d2)
        holo__domain_interface_area = interdomain_surface_computer(s2d1, s2d2)

        # *2 = 218, 933 vs paper -- 288, 1024  # oni to asi nedělej dvěma..., ale priblibzne to odpovida
        self.assertAlmostEqual(288, 2*apo__domain_interface_area, delta=0.3 * 288)  # 218
        self.assertAlmostEqual(1024, 2*holo__domain_interface_area, delta=0.3 * 1024)  # 933

        # rmsd
        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

        print(get_rmsd(s1d1+s1d2, s2d1+s2d2))  # 10.1 vs 8.0 in paper, todo celkem velky rozdil...

        # hinge
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
        angle, translation = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)
        self.assertAlmostEqual(147, 180/np.pi * angle, delta=2)

    def test_get_sasa_for_structure(self):
        s = self.get_test_structure()
        chain_b = s[0]['B']
        residues = ChainResidues(list(chain_b), s.id, 'B')

        sasa_computer = GetSASAForStructure()
        chain_b_sasa = sasa_computer(residues)

        self.assertGreater(chain_b_sasa, 1)

    def test_rmsd_guanylate_kinase_paper(self):
        apo = self.load_test_structure('1ex6')
        holo = self.load_test_structure('1ex7')

        apo_chain = apo[0]['B']  # note that different chain (as by dyndom), why?
        holo_chain = holo[0]['A']

        logging.root.setLevel(logging.INFO)
        self.assertTrue(sequences_same(apo_chain, holo_chain))

        # analyze just A chain (like in the paper)
        apo_residues = ChainResidues.from_bio_chain(apo_chain)
        holo_residues = ChainResidues.from_bio_chain(holo_chain)


        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rmsd = GetRMSD((get_centered_c_alpha_coords, GetRotationMatrix((get_centered_c_alpha_coords,))))

        rmsd = get_rmsd(apo_residues, holo_residues)
        self.assertAlmostEqual(4.4, rmsd, delta=0.1)  # 4.37 (vs 4.4 Å in the paper)

    def test_hinge_guanylate_kinase_paper(self):
        apo = self.load_test_structure('1ex6')
        holo = self.load_test_structure('1ex7')

        apo_chain = apo[0]['B']  # note that different chain (as by dyndom), why?
        holo_chain = holo[0]['A']

        logging.root.setLevel(logging.INFO)
        self.assertTrue(sequences_same(apo_chain, holo_chain))

        apo_d1 = DomainResidues.from_domain(DomainResidueMapping('D1', 'B', [200 + 1, 200 + 84], [200 + 32, 200 + 186]), apo)
        apo_d2 = DomainResidues.from_domain(DomainResidueMapping('D2', 'B', [200 + 33], [200 + 83]), apo)

        # zmena chainu preci nepomahala, tak kde je zakopany pes?

        holo_d1 = DomainResidues.from_domain(DomainResidueMapping('D1', 'A', [1, 84], [32, 186]), holo)
        holo_d2 = DomainResidues.from_domain(DomainResidueMapping('D2', 'A', [33], [83]), holo)

        get_c_alpha_coords = GetCAlphaCoords()
        get_centroid = GetCentroid((get_c_alpha_coords,))
        get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
        get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))
        get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))

        angle, translation_in_axis = get_hinge_angle(apo_d1, apo_d2, holo_d1, holo_d2)
        self.assertAlmostEqual(47, 180/np.pi*angle, delta=0.2)  # in paper: dyndom: 47°, their principal axes 43.9

    def test_rmsd_pheromone_binding_protein_paper(self):
        apo = self.load_test_structure('2fjy')
        holo = self.load_test_structure('1dqe')

        apo_chain = apo[0]['A']
        holo_chain = holo[0]['A']

        # # logging.root.setLevel(logging.INFO)
        # self.assertTrue(sequences_same(apo_chain, holo_chain))

        # chains also have leading and trailing residues not present in the other, remove them

        apo_pp = chain_to_polypeptide(apo_chain)[:-5]
        holo_pp = chain_to_polypeptide(holo_chain)[6:]

        # analyze just A chain (like in the paper)
        apo_residues = ChainResidues(list(apo_pp), apo.id, apo_chain.id)
        holo_residues = ChainResidues(list(holo_pp), holo.id, holo_chain.id)

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

        interdomain_surface_computer = GetInterdomainSurface((GetSASAForStructure(),))
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

        # test rmsd is still 0
        chain_a = ChainResidues.from_bio_chain(s1[0]['A'])
        chain_a_rotated = ChainResidues.from_bio_chain(s2[0]['A'])

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

        s1d1 = ChainResidues.from_bio_chain(s1[0]['A'])
        s1d2 = ChainResidues.from_bio_chain(s1[0]['B'])
        s2d1 = ChainResidues.from_bio_chain(s2[0]['A'])
        s2d2 = ChainResidues.from_bio_chain(s2[0]['B'])

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

        angle, translation_in_axis = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)

        self.assertAlmostEqual(ANGLE, angle, places=3)
        self.assertAlmostEqual(TRANSLATION_IN_AXIS, translation_in_axis, places=3)
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

        s1d1 = ChainResidues.from_bio_chain(s1[0]['A'])
        s1d2 = ChainResidues.from_bio_chain(s1[0]['B'])
        s2d1 = ChainResidues.from_bio_chain(s2[0]['A'])
        s2d2 = ChainResidues.from_bio_chain(s2[0]['B'])

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

        angle, translation_in_axis = get_hinge_angle(s1d1, s1d2, s2d1, s2d2)

        self.assertAlmostEqual(ANGLE, angle, places=3)
        self.assertAlmostEqual(TRANSLATION_IN_AXIS, translation_in_axis, places=3)
        # GetHingeAngle does not return AXIS_DIRECTION, or AXIS_LOCATION (yet)
