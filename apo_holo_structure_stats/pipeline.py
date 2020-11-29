from Bio.PDB import MMCIFParser

from apo_holo_structure_stats.input.download import get_structure_stringio


def single_struct_analyzer(struct):
    pass


def parse_stucture(structure_code):
    structure_file = get_structure_stringio(structure_code)

    mmcif_parser = MMCIFParser()

    structure = mmcif_parser.get_structure(structure_code, structure_file)
    mmcif_dict_undocumented = mmcif_parser._mmcif_dict

    return structure


def main():
    structure_group = []
    apo, holo = filter(, structure_group)

    # gets mapping between residues
    apo_holo_pairs = make_pairs
    holo_holo_pairs = make_pairs

    for pair in apo_holo_pairs:
        twinstructanalyzer(pair)

    for pair in holo_holo_pairs:
        twinstructanalyzer(pair)
