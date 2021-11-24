""" List all structures from apo_holo.dat as single line of csv, so it can be pasted as an argument to
pipeline/filter_structures.py script; for testing purposes (paper repl).
"""

import pandas as pd

from apo_holo_structure_stats.paper_repl.main import get_paper_apo_holo_dataframe


if __name__ == '__main__':
    df = get_paper_apo_holo_dataframe()

    structures = []
    for index, row in df.iterrows():
        s1_pdb_code = row.apo[:4]
        s2_pdb_code = row.holo[:4]
        structures.append(s1_pdb_code)
        structures.append(s2_pdb_code)

    print(len(structures))
    print(repr(structures))
    print(','.join(structures))
