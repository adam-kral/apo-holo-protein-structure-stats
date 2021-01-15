from pathlib import Path
from unittest import TestCase

from Bio.PDB import MMCIFParser

from apo_holo_structure_stats.core.analyses import SetOfResidues, GetSASAForStructure, GetInterdomainSurface, ChainResidues


class TestAnalyses(TestCase):
    def setUp(self) -> None:
        self.structure = MMCIFParser().get_structure('2gdu', Path(__file__).parent / 'test_data' / '2gdu.cif')

    def get_test_structure(self):
        return self.structure.copy()

    def test_get_sasa_for_structure(self):
        s = self.get_test_structure()
        chain_b = s[0]['B']
        residues = ChainResidues(list(chain_b), s.id, 'B')

        sasa_computer = GetSASAForStructure()
        chain_b_sasa = sasa_computer(residues)
        print(chain_b_sasa)

        self.assertGreater(chain_b_sasa, 1)

    def test_interdomain_surface(self):
        s = self.get_test_structure()
        chain_a = ChainResidues(list(s[0]['A']), s.id, 'A')
        chain_b = ChainResidues(list(s[0]['B']), s.id, 'B')

        interdomain_surface_computer = GetInterdomainSurface((GetSASAForStructure(),))
        area = interdomain_surface_computer(chain_a, chain_b)
        print(area)
        self.assertGreater(area, 1)
