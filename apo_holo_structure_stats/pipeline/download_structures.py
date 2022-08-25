import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Iterable

import pandas as pd

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.input.download import find_or_download_structure
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args, get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks

logger = logging.getLogger(__name__)

MAX_FUTURES = 2000  # max Future objects at once (take up memory)


def download_structures(pdb_codes: List[str], workers: int) -> Iterable[Path]:
    with ThreadPoolExecutor(workers) as executor:
        futures = submit_tasks(executor, MAX_FUTURES, find_or_download_structure, pdb_codes)
        for i, (pdb_code, future) in enumerate(zip(pdb_codes, futures)):
            structure_path = future.result()
            logger.info(f'Downloaded structure {pdb_code}, {i+1}/{len(pdb_codes)} total')


def main():
    import argparse
    import sys

    parser = get_argument_parser()
    parser.add_argument('--download_threads', type=int, default=1, help='number of threads')

    # parser.add_argument('--pdb_dir', type=str, action='store_true',
    #                     help='now pdb_codes_or_directory is a path to a directory with mmcif files. Whole tree structure is inspected, all files are assumed to be mmcifs.')
    parser.add_argument('-i', '--input_type', default='json', choices=['json', 'pdb_codes'], help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('input', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
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

    download_structures(list(pdb_codes), args.download_threads)


if __name__ == '__main__':
    main()
