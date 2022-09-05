""" Download structures from the PDB.

See `ah-download-structures --help` for details and usage.
"""
import argparse
import itertools
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Iterable

import more_itertools
import pandas as pd

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.input.download import find_or_download_structure
from apo_holo_structure_stats.pipeline.utils.log import get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks

logger = logging.getLogger(__name__)

MAX_FUTURES = 2000  # max Future objects at once (take up memory)


def download_structures(pdb_codes: List[str], workers: int) -> Iterable[Path]:
    """ Download structures from the PDB.

    Note it returns a generator of paths to the downloaded structures which __must be
    consumed__ to actually download the structures. todo return a list instead?"""
    with ThreadPoolExecutor(workers) as executor:
        futures = submit_tasks(executor, MAX_FUTURES, find_or_download_structure, pdb_codes)
        for i, (pdb_code, future) in enumerate(zip(pdb_codes, futures)):
            structure_path = future.result()
            yield structure_path
            logger.info(f'Downloaded structure {pdb_code}, {i+1}/{len(pdb_codes)} total')


parser = get_argument_parser(description='Download structures from the PDB.',
                             formatter_class=argparse.RawDescriptionHelpFormatter,
                             epilog='''\
Files will be downloaded to the Settings.STRUCTURE_STORAGE_DIRECTORY.
Other scripts will automatically use this directory for loading the structures.

Usage:
    ah-download-structures -v --workers 10 chains.json  
     ah-download-structures -v -i pdb_codes 1abc,2abc                                             
''')
parser.add_argument('--threads', type=int, default=1, help='Number of download threads.')
parser.add_argument('-i', '--input_type', default='json', choices=['json', 'pdb_codes'], help='Set `pdb_codes` if the input is a csv list of PDB codes.')
parser.add_argument('input', help='JSON (records) dataset of chains, required columns: pdb_code, chain_id. '
                                  'Or if --input_type pdb_codes option is present, a csv string of PDB codes.')


def main():
    import sys
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()

    if args.input_type == 'pdb_codes':
        pdb_codes = args.input.strip().split(',')

        if not pdb_codes:
            logger.error('No pdb codes specified')
            sys.exit(1)

    elif args.input_type == 'json':
        # todo to accept just list of pdb_codes, add column chain_id? And won't that break chain_whitelists (want empty set)
        chains = pd.read_json(args.input)
        chains_gb_pdb_code = chains.groupby('pdb_code')
        pdb_codes = chains_gb_pdb_code.indices.keys()
    else:
        raise ValueError(f'Unknown input_type: {args.input_type}')

    logger.info(f'Downloading {len(pdb_codes)} structures.')
    more_itertools.consume(
        download_structures(list(pdb_codes), args.threads)
    )


if __name__ == '__main__':
    main()
