import logging

import pandas as pd
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model

from apo_holo_structure_stats.core.analysesinstances import get_chains


def get_paper_apo_holo_dataframe(filename='apo_holo.dat'):
    df = pd.read_csv(filename, delimiter=r'\s+', comment='#', header=None,
                     names=('apo_id', 'holo_id', 'domain_count', 'ligand_codes'), dtype={'domain_count': int})

    # split pdb_code and chain id (chain id can be _)
    df['apo_chain_id'] = df['apo_id'].map(lambda s: s[4:])
    df['holo_chain_id'] = df['holo_id'].map(lambda s: s[4:])
    df['apo_pdb_code'] = df['apo_id'].map(lambda s: s[:4])
    df['holo_pdb_code'] = df['holo_id'].map(lambda s: s[:4])
    # drop now superfluous columns
    df = df.drop(columns=['apo_id', 'holo_id'])
    return df


def get_chain_by_chain_code(model: Model, paper_chain_code: str) -> Chain:
    if paper_chain_code == '_':
        # seems that this is the code when there's only a single chain (probably named A?)
        chain = next(iter(model))

        if len(model) != 1 or chain.id != 'A':
            logging.warning(f'{model.get_parent().id} {paper_chain_code}, underscore is not what I originally thought. {chain.id}'
                            f', {len(model)}')
            long_chains = get_chains(model)
            if len(long_chains) > 1:
                logging.warning(f'model contains {len(long_chains)} chains with 50+ aa, underscore is not what I think it is now')
            return long_chains[0]  # sometimes there is also chain B with ligands.

        return chain
    return model[paper_chain_code]