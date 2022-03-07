import argparse
import logging
from pathlib import Path
from typing import Set

import pandas as pd

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.input.uniprot_groups import get_basic_uniprot_groups__df
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args

logger = logging.getLogger(__name__)

OUTPUT_DIR = 'output_dir'  # hack, move to settings, or shell args


def collect_chains_for_uniprot_ids(uniprot_ids: Set[str] = None, limit_size: int = None, seed: int = None):
    logger.info('Downloading and parsing a csv listing all chains in PDB...')

    Path(OUTPUT_DIR).mkdir(exist_ok=True)
    chains = get_basic_uniprot_groups__df(download_dir=OUTPUT_DIR)
    # groups = chains.groupby('uniprotkb_id').size()

    if uniprot_ids:
        chains = chains[chains.uniprotkb_id.isin(uniprot_ids)]

    if limit_size:
        ok_chains = chains[chains.uniprot_group_size <= limit_size]
        chains_from_oversized_groups = chains[chains.uniprot_group_size > limit_size]

        # subsample chains_from_oversized_groups
        assert seed is not None
        subsampled = chains_from_oversized_groups.groupby('uniprotkb_id').sample(limit_size, random_state=seed)

        chains = pd.concat([ok_chains, subsampled])

    return chains


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--uniprot_ids', help='comma-separated list of primary uniprot accessions')
    parser.add_argument('--chains', help='json (records) dataset of chains that will be augmented with uniprotkb_id')
    parser.add_argument('--limit_group_size_to', type=int, help='comma-separated list of primary uniprot accessions')
    parser.add_argument('--seed', default=42, type=int, help='comma-separated list of primary uniprot accessions')
    parser.add_argument('output_file', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')
    add_loglevel_args(parser)

    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    if args.uniprot_ids:
        uniprot_ids = set((g.strip() for g in args.uniprot_ids.split(',')))
    else:
        # means all uniprot groups
        uniprot_ids = None

    chains = collect_chains_for_uniprot_ids(uniprot_ids, args.limit_group_size_to, args.seed)

    if args.chains:
        # chain dataset was input, only add uniprotkb_id field
        input_chains = pd.read_json(args.chains)
        chains = input_chains.merge(chains, how='left', on=['pdb_code', 'chain_id'])

        chains_without_uniprot = chains[pd.isna(chains.uniprotkb_id)]
        if not chains_without_uniprot.empty:
            logger.warning('UniprotKB ID not found for some chains, won`t be included in the result:\n'
                           + str(chains_without_uniprot))
            chains = chains.dropna(subset=['uniprotkb_id'])

    chains.to_json(args.output_file, orient='records')


if __name__ == '__main__':
    main()
