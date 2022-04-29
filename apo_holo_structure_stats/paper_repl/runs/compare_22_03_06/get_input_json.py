import logging
from pathlib import Path

import pandas as pd

from apo_holo_structure_stats import settings
from apo_holo_structure_stats.input.download import parse_mmcif
from apo_holo_structure_stats.paper_repl.main import get_paper_apo_holo_dataframe, get_chain_by_chain_code

df = get_paper_apo_holo_dataframe('/home/adam/pschool/bakalarka/apo-holo-protein-structure-stats/manual_test_data/apo_holo.dat')

apo_chains = df[['apo_pdb_code', 'apo_chain_id']].rename(columns={'apo_pdb_code': 'pdb_code', 'apo_chain_id': 'chain_id'})
holo_chains = df[['holo_pdb_code', 'holo_chain_id']].rename(columns={'holo_pdb_code': 'pdb_code', 'holo_chain_id': 'chain_id'})
chains = pd.concat([apo_chains, holo_chains])

# fix the underscores in the original dataset:

def get_actual_chain_id(s1_pdb_code: str, s1_paper_chain_code: str):
    logging.info(f'processing {s1_pdb_code} {s1_paper_chain_code}')
    apo_parsed = parse_mmcif(s1_pdb_code)

    apo = apo_parsed.structure
    apo_residue_id_mappings = apo_parsed.bio_to_mmcif_mappings
    # apo_poly_seqs = apo_parsed.poly_seqs

    # get the first model (s[0]), (in x-ray structures there is probably a single model)
    apo, apo_residue_id_mappings = map(lambda s: s[0], (apo, apo_residue_id_mappings))

    apo_chain = get_chain_by_chain_code(apo, s1_paper_chain_code)
    actual_chain_id = apo_chain.id  # not an underscore (perhaps standing for implicit chain id)
    return actual_chain_id

# logging.basicConfig()
logging.root.setLevel(logging.INFO)

settings.STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path('../../pdb_structs')
chains['chain_id'] = list(map(lambda row: get_actual_chain_id(row.pdb_code, row.chain_id), chains.itertuples()))
print(chains)
chains.to_json('chains.json', orient='records')



