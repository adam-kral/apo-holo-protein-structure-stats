#!/usr/bin/env python3
import json

from apo_holo_structure_stats.input.download import get_best_isoform_for_chains


def get_isoform(pdb_code: str, chain_id: str) -> str:
    """ Returns uniprot accession for the best isoform mapped to the structure's chain. """

    isoform_mapping = get_best_isoform_for_chains(pdb_code)
    return isoform_mapping[chain_id][0].uniprot_id


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('json_file', help='annotate the list of structures with isoform data. File needs to contain list of objects with pdb_code and main_chain_id')
    args = parser.parse_args()

    with open(args.json_file) as f:
        structures_info = json.load(f)

    for s in structures_info:
        s['isoform_id'] = get_isoform(s['pdb_code'], s['main_chain_id'])

    with open(args.json_file, 'w') as f:
        json.dump(structures_info, f)
