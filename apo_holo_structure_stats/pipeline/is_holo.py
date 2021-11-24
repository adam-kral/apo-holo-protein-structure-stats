#!/usr/bin/env python3
import concurrent.futures
import itertools
import json
import logging
import warnings
from typing import Dict

from Bio.PDB import MMCIFParser

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.pipeline.log import add_loglevel_args
from apo_holo_structure_stats.core.analyses import IsHolo

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=1, help='number of subprocesses')
    parser.add_argument('structures_json', help='annotate the list of structures with is_holo bool. File needs to contain list of objects with pdb_code and chain_id and path to the structure')
    parser.add_argument('output_file', help='writes input json annotated with boolean "is holo"')
    add_loglevel_args(parser)
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    # todo multiprocessing logging will only work to stdout: "Although logging is thread-safe, and logging to a single file from multiple threads in a single process is supported, logging to a single file from multiple processes is not supported, because there is no standard way to serialize access to a single file across multiple processes in Python.
    # todo ma-li struktura vic chainu, nacist jen jednou (ale budu to vsechno delat po skupinach, tak preci musim vedet, jaky chain tak nejak chci, tak by byl jen ten. Ale muze jich byt vic..)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        def is_holo(ordinal, s: Dict):
            logger.info(f'processing {ordinal}-th structure {s["pdb_code"]}')

            is_holo_analyzer = IsHolo()

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = MMCIFParser().get_structure(s['pdb_code'], s['path'])[0]

            return is_holo_analyzer(model, model[s['chain_id']])

        results = executor.map(is_holo, itertools.count(), structures_info)

    # add is_holo flag to the structure info dicts
    for is_holo, structure_info in zip(results, structures_info):
        structure_info['is_holo'] = is_holo

    with open(args.output_file, 'w') as f:
        json.dump(structures_info, f)
