""" Obtain domains and secondary structure for (the already paired) structures.

This script obtains it, given the pdb_codes in the pairs json,  using the pdbe-kb API.
It is not extensible (but could be), currently
users are expected to use their data gathering scripts to obtain additional data they need in run_analyses.py.

There are much fewer apo-holo paired structures than in the whole pdb and the APIs "rely on user restraint".
"""

import itertools
import logging
import shelve
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analyses import GetSecondaryStructureForStructure, GetDomainsForStructure
from apo_holo_structure_stats.pipeline.make_pairs_lcs import load_pairs_json, pairs_without_mismatches
from apo_holo_structure_stats.pipeline.utils.json import maybe_print
from apo_holo_structure_stats.pipeline.utils.log import get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks

logger = logging.getLogger(__name__)

get_ss = GetSecondaryStructureForStructure()
get_domains = GetDomainsForStructure()


parser = get_argument_parser(description=__doc__)
parser.add_argument('--threads', default=12, type=int, help='API download threads.')
parser.add_argument('pairs_json', help='Output of ah-make-pairs.')


def main():
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()
    # the pairs json could be already pairs without mismatches, I would save in the script 2/3 of the pairs, memory..

    potential_pairs = load_pairs_json(args.pairs_json)
    if potential_pairs.empty:
        logger.warning('Input json contains no records.')
        sys.exit(0)

    print(potential_pairs)
    pairs = pairs_without_mismatches(potential_pairs)
    # pairs = pairs.iloc[:40]  # todo testing hack
    print(pairs)

    def rename_col(old_name):
        cols = ['pdb_code', 'chain_id']
        for col in cols:
            if old_name.startswith(col):
                return col
        raise ValueError('column not found')

    all_chains = pd.concat([
        pairs[['pdb_code_apo', 'chain_id_apo']].rename(columns=rename_col),
        pairs[['pdb_code_holo', 'chain_id_holo']].rename(columns=rename_col),
    ])

    all_structs = all_chains['pdb_code'].unique()
    print(all_structs)

    quiet = args.loglevel > logging.INFO

    # could open the shelf with flag fast and file sync once in a while, but API requests are slow anyway
    with ThreadPoolExecutor(max_workers=args.threads) as executor:

        logger.info(f'total structs: {len(all_structs)}')
        successes = errors = 0

        ss_futures = submit_tasks(executor, 40 * args.threads, get_ss, all_structs)

        with shelve.open('db_get_ss') as db:
            for i, pdb_code, f in zip(itertools.count(), all_structs, ss_futures):
                try:
                    db[pdb_code] = f.result()
                    successes += 1
                except Exception:
                    logger.exception(f'get_ss for {pdb_code} failed with:')
                    errors += 1

                maybe_print(quiet, f'\r suc {successes}, err {errors} done {i+1}/{len(all_structs)}', end='')

        successes = errors = 0

        domain_futures = submit_tasks(executor, 40 * args.threads, get_domains, all_structs)

        with shelve.open('db_get_domains') as db:
            for i, pdb_code, f in zip(itertools.count(), all_structs, domain_futures):
                try:
                    db[pdb_code] = f.result()
                    successes += 1
                except Exception:
                    logger.exception(f'get_domains for {pdb_code} failed with:')
                    errors += 1

                maybe_print(quiet, f'\r suc {successes}, err {errors} done {i+1}/{len(all_structs)}', end='')


if __name__ == '__main__':
    main()
