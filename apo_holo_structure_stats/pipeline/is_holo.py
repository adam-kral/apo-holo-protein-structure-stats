#!/usr/bin/env python3
import concurrent.futures
import itertools
import json
import logging
import warnings

from Bio.PDB import MMCIFParser

from apo_holo_structure_stats.pipeline.log import add_loglevel_args
from apo_holo_structure_stats.core.analyses import IsHolo, GetMainChain, GetChains


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=1, help='number of subprocesses')
    parser.add_argument('structures_json', help='annotate the list of structures with is_holo bool. File needs to contain list of objects with pdb_code and main_chain_id and path to the structure')
    parser.add_argument('output_file', help='writes input json annotated with boolean "is holo"')
    add_loglevel_args(parser)
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        def is_holo(ordinal, s):
            logging.info(f'processing {ordinal}-th structure {s["pdb_code"]}')

            is_holo_analyzer = IsHolo((lambda model: model[s['main_chain_id']],))  # tady jsem dal do dependecies lambdu místo GetMainChain (taky tim porušuju typ argumentu)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = MMCIFParser().get_structure(s['pdb_code'], s['path'])[0]

            return is_holo_analyzer(model)

        results = executor.map(is_holo, itertools.count(), structures_info)

    # add is_holo flag to the structure info dicts
    for is_holo, structure_info in zip(results, structures_info):
        structure_info['is_holo'] = is_holo

    with open(args.output_file, 'w') as f:
        json.dump(structures_info, f)
