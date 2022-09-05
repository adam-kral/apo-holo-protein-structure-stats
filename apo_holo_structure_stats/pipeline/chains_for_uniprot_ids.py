""" Collect PDB chains with their uniprot ids.

This script downloads a csv file listing all chains in PDB and returns a json file with chains with uniprot ids.
See `ah-chains-uniprot --help` for details.
"""
import argparse
import logging
import sys
from pathlib import Path
from typing import Set

import pandas as pd

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.input.uniprot_groups import get_chains_with_uniprot_ids
from apo_holo_structure_stats.pipeline.utils.log import get_argument_parser

logger = logging.getLogger(__name__)


def collect_chains_for_uniprot_ids(uniprot_ids: Set[str] = None, limit_size: int = None, seed: int = None):
    logger.info('Downloading and parsing a csv listing all chains in PDB...')

    chains = get_chains_with_uniprot_ids(download_dir=Path())

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


parser = get_argument_parser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description='Collect PDB chains with their uniprot ids.',
                             epilog='''\
By default all PDB chains are collected (which are in the SIFTS service).
Output fields are: pdb_code, chain_id, uniprotkb_id, uniprot_group_size
    where uniprot_group_size is the number of chains in the PDB that have the same uniprot id.

Data are obtained from SIFTS' uniprot_segments_observed.csv file.

Usage:
    ah-chains-uniprot chains.json
    ah-chains-uniprot --chains <chains_without_uniprot>.json chains.json
    ah-chains-uniprot --uniprot_ids P12345,P12346 chains.json
    ah-chains-uniprot --limit_group_size_to 10 --seed 42 chains.json
''')
parser.add_argument('--uniprot_ids', help='If set, output only chains with these uniprots\n'
                                          'comma-separated list of primary uniprot accessions')
parser.add_argument('--chains', help='If set, this json (records) dataset of chains will be augmented '
                                     'with uniprotkb_id and uniprot_group_size column. Required columns: '
                                     'pdb_code, chain_id')

parser.add_argument('--limit_group_size_to', type=int, help='If set, all chains from uniprot groups of this size or '
                                                            'smaller will be included. Chains from larger groups '
                                                            'will be subsampled using --seed')
parser.add_argument('--seed', default=42, type=int, help='Seed for subsampling chains. See --limit_group_size_to')

parser.add_argument('output_file', help='Output file path. JSON (records) dataset of chains with columns: '
                                        'pdb_code, chain_id, uniprotkb_id, uniprot_group_size')


def main():
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()

    if args.uniprot_ids:
        uniprot_ids = set((g.strip() for g in args.uniprot_ids.split(',')))
    else:
        # means all uniprot groups
        uniprot_ids = None

    chains = collect_chains_for_uniprot_ids(uniprot_ids, args.limit_group_size_to, args.seed)

    if args.chains:
        # chain dataset was input, only add uniprotkb_id field
        try:
            input_chains = pd.read_json(args.chains)
        except ValueError:
            sys.stderr.write('Error: input file is not a valid json dataset of chains.')
            sys.exit(1)

        assert 'pdb_code' in input_chains.columns and 'chain_id' in input_chains.columns

        chains = input_chains.merge(chains, how='left', on=['pdb_code', 'chain_id'])

        chains_without_uniprot = chains[pd.isna(chains.uniprotkb_id)]
        if not chains_without_uniprot.empty:
            logger.warning('UniprotKB ID not found for some chains, won`t be included in the result:\n'
                           + str(chains_without_uniprot))
            chains = chains.dropna(subset=['uniprotkb_id'])

    chains.to_json(args.output_file, orient='records')


if __name__ == '__main__':
    main()
