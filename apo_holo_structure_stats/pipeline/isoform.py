#!/usr/bin/env python3
import concurrent.futures
import itertools
import json
import logging

from apo_holo_structure_stats.input.download import get_best_isoform_for_chains, APIException
from apo_holo_structure_stats.pipeline.log import add_loglevel_args


def get_isoform(pdb_code: str, chain_id: str) -> str:
    """ Returns uniprot accession for the best isoform mapped to the structure's chain. """

    isoform_mapping = get_best_isoform_for_chains(pdb_code)
    return isoform_mapping[chain_id][0].uniprot_id


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=1, help='number of threads for concurrent API requests')
    parser.add_argument('structures_json', help='annotate the list of structures with isoform data. File needs to contain list of objects with pdb_code and chain_id')
    parser.add_argument('output_file', help='writes input json annotated with isoform uniprot id')
    add_loglevel_args(parser)
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        def get_isoform_or_none(ordinal, s):
            logging.info(f'processing {ordinal}-th structure {s["pdb_code"]}')

            try:
                return get_isoform(s['pdb_code'], s['chain_id'])
            except APIException as e:
                if '404' in str(e.__cause__):  # todo thihs does not work EDIT it does catch 404, so why I wrote that?
                    logging.info(f'isoform not found for {s["pdb_code"]}')
                else:
                    logging.exception(f'api error for {s["pdb_code"]}')
            except Exception:
                logging.exception(f'unexpected error for {s["pdb_code"]}')

            return None

        # todo if more chains, for a structure, cache (problem with multithreading though)
        results = executor.map(get_isoform_or_none, itertools.count(), structures_info)

    # update the dict (json) with isoform information
    for isoform, structure_info in zip(results, structures_info):
        structure_info['isoform'] = isoform

    with open(args.output_file, 'w') as f:
        json.dump(structures_info, f)
