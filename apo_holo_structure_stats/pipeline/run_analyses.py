#!/usr/bin/env python3
import itertools
import json
from typing import List


def run_analyses_for_isoform_group(apo_codes: List[str], holo_codes: List[str], get_structure):
    # apo-holo analyses
    apo_holo_pairs = itertools.product(apo_codes, holo_codes)

    for apo_code, holo_code in apo_holo_pairs:
        print(apo.id, holo.id)
        print(rmsd_if_seqs_same(apo, holo))

    # holo-holo analyses
    holo_holo_pairs = itertools.combinations(holo_codes, 2)

    for holo1_code, holo2_code in holo_holo_pairs:
        print(holo1.id, holo2.id)
        print(rmsd_if_seqs_same(holo1, holo2))


if __name__ == '__main__':
    # run for all isoforms default
    # optionally specify a single isoform

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--isoform', help='process only structures with main chain of that isoform')
    parser.add_argument('structures_json_file', help='list of structures {pdb_code: , path: , isoform_id: , is_holo: bool, ?main_chain_id: }')
    parser.add_argument('output_file', help='dumped results of analyses')
    args = parser.parse_args()

    with open(args.structures_json_file) as f:
        structures_info = json.load(f)

    if args.isoform is not None:
        structures_info = list(filter(lambda s: s['isoform_id'] == args.isoform, structures_info))


    with open(args.output_file, 'w') as f:
        # serializace analyz bude asi prubezna

        # groupby isoform
            # run_analyses_for_isoform_group(, get_structure = lambda: 1  )
        pass
