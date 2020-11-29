import unittest

import pandas as pd

from ..input.uniprot_groups import hash_dataframe_rows_no_categorize


class TestUniprotGroups(unittest.TestCase):
    def test_string_hashing_returns_same_hash_even_if_string_are_different_in_memory(self):
        """ pandas.read_csv apparently internalizes strings; so it could be deceiving that it works, but sometimes, when
        two identical strings wouldn't point on an identical object in memory (pointers in pandas-underlying numpy array),
        it could stop working
        """

        s1 = 'hovno'
        s2 = s1[:2] + s1[2:]

        self.assertIsNot(s1, s2)

        df = pd.DataFrame({'col1': [1, 1], 'col2': [s1, s2]})

        row_hashes = hash_dataframe_rows_no_categorize(df)
        self.assertEqual(row_hashes[0], row_hashes[1])


if __name__ == '__main__':
    unittest.main()
