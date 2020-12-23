#!/usr/bin/env python3
import json

from Bio.PDB import MMCIFParser

from apo_holo_structure_stats.core.analyses import IsHolo, GetMainChain, GetChains


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('structures_json', help='annotate the list of structures with is_holo bool. File needs to contain list of objects with pdb_code and main_chain_id and path to the structure')
    parser.add_argument('output_file', help='writes input json annotated with boolean "is holo"')
    args = parser.parse_args()

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    # add is_holo flag to the structure info dicts
    for s in structures_info:
        is_holo_analyzer = IsHolo((lambda model: model[s['main_chain_id']],))  # tady jsem dal do dependecies lambdu místo GetMainChain (taky tim porušuju typ argumentu). Ale udelat to zde muzu, kdyz uz mam main_chain_id, todo mozna napsat do helpu json_file argu, ze tam musi byt i main_chain_id

        model =  MMCIFParser().get_structure(s['pdb_code'], s['path'])[0]

        s['is_holo'] = is_holo_analyzer(model)

    with open(args.output_file, 'w') as f:
        json.dump(structures_info, f)
