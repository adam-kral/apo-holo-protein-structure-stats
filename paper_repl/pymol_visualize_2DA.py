import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analysesinstances import get_rotation_matrix, get_centroid
from apo_holo_structure_stats.core.dataclasses import SetOfResidues, DomainResidueMapping, DomainResidues
from apo_holo_structure_stats.core.json_deserialize import unfold_tuple_to_columns, tuple_it
from apo_holo_structure_stats.input.download import parse_mmcif

from .main import get_longest_common_polypeptide  # todo use the one from apo_holo_structure_stats, main is
# deprecated and can be deleted


OUTPUT_DIR = 'output'
logger = logging.getLogger(__name__)


def superimpose_structure(s1: Model, by_1: SetOfResidues, on_2: SetOfResidues):
    """ Rotate and move the whole s1 structure on a second structure. So that by_1 will be superimposed onto on_2.

    by_1 and on_2 can be a subset of the structure's residues (e.g. domains). Both have to have the same length and
    should contain corresponding (sequence aligned) residues, in the same order.
    :returns a copy of s1 Model with all its atom coordinates rotated and translated (by the superposition).
    """
    s1_copy = s1.copy()

    # put all coordinates of s1 into a np.array
    # (get_unpacked_list [in contrary to Model.get_atoms()] also includes disordered atoms/residues. So this would
    # handle that correctly, however we skip structures with disordered residues anyway in our analyses.)
    s1_atoms = [a for c in s1_copy.get_chains() for r in c.get_unpacked_list() for a in r.get_unpacked_list()]
    s1_coords = np.array([atom.get_coord() for atom in s1_atoms])

    # superimpose the two structures by their first domain
    U = get_rotation_matrix(by_1, on_2)
    new_s1_coords = (s1_coords - get_centroid(by_1)) @ U + get_centroid(on_2)

    for atom, new_coord in zip(s1_atoms, new_s1_coords):
        atom.set_coord(new_coord)

    return s1_copy


@dataclass
class TwoDomainArrangement:
    pdb_code: str
    d1: DomainResidueMapping
    d2: DomainResidueMapping


def visualize_2DA(apo_2DA, holo_2DA, paper_apo_spans):
    """ Writes superimposed holo structure to a file, prints Pymol script which can be directly pasted in pymol.

     Printed Pymol script will:
     1) automatically load both structures (superimposed holo from filesystem, apo from the internet)
     2) create objects and selections for domains, and the two-domain arrangements
     3) color the selections by domain, apo/holo and paper/ours
        - colors - ours more saturation, paper faded
            - red, yellow apo (first and second domain respectively)
            - green, blue holo
     4) provide example usage in the last script paragraph
     """

    # load the structure from file
    a = parse_mmcif(apo_2DA.pdb_code)
    h = parse_mmcif(holo_2DA.pdb_code)
    apo = a.structure
    holo = h.structure

    ###### vlozene z mainu
    apo_mapping = a.bio_to_mmcif_mappings[0][apo_2DA.d1.chain_id]
    holo_mapping = h.bio_to_mmcif_mappings[0][holo_2DA.d1.chain_id]

    # crop polypeptides to longest common substring
    c1_common_seq, c2_common_seq = get_longest_common_polypeptide(a.poly_seqs[apo_mapping.entity_poly_id], h.poly_seqs[holo_mapping.entity_poly_id])
    c1_label_seq_ids = list(c1_common_seq.keys())
    c2_label_seq_ids = list(c2_common_seq.keys())

    label_seq_id_offset = c2_label_seq_ids[0] - c1_label_seq_ids[0]
    ###### end vlozene

    # get residues of the first domain, in both apo and holo structures
    apo_d1 = DomainResidues.from_domain(apo_2DA.d1, apo[0], apo_mapping)
    holo_d1 = DomainResidues.from_domain(holo_2DA.d1, holo[0], holo_mapping)
    # superimpose holo onto apo, using the first domain
    superimposed_holo_model = superimpose_structure(holo[0], holo_d1, apo_d1)
    # save the structure
    name = holo.id + f'_{holo_d1.domain_id}onto_{apo_d1.domain_id}'
    io = MMCIFIO()
    superimposed_holo = Structure(name)
    superimposed_holo.add(superimposed_holo_model)
    io.set_structure(superimposed_holo)
    sholo_file_path = Path(OUTPUT_DIR, name + '.cif')
    io.save(str(sholo_file_path), preserve_atom_numbering=True)

    def get_resi_selection(spans):
        selection = []
        for from_, to in spans:
            selection.append(f'resi {from_}-{to}')

        return '(' + ' or '.join(selection) + ')'

    # convert paper spans to label seqs, so we can show them in Pymol
    def get_paper_domain(d: DomainResidueMapping, paper_spans, residue_id_mapping):
        # translate spans to label seq ids and return a domain object
        segment_beginnings = list(map(residue_id_mapping.find_label_seq, np.array(paper_spans)[:, 0].tolist()))
        segment_ends = list(map(residue_id_mapping.find_label_seq, np.array(paper_spans)[:, 1].tolist()))
        logger.debug(segment_beginnings)
        logger.debug(segment_ends)
        return DomainResidueMapping(d.domain_id, d.chain_id, segment_beginnings, segment_ends)

    logger.debug(paper_apo_spans)  # [d1, d2] where d1 [(), (),...]
    paper_apo_drm1 = get_paper_domain(apo_2DA.d1, paper_apo_spans[0], apo_mapping)
    paper_apo_drm2 = get_paper_domain(apo_2DA.d2, paper_apo_spans[1], apo_mapping)
    label_seq_id_offset = c2_label_seq_ids[0] - c1_label_seq_ids[0]
    paper_holo_drm1 = DomainResidueMapping.from_domain_on_another_chain(paper_apo_drm1, holo_d1.chain_id, label_seq_id_offset)
    paper_holo_drm2 = DomainResidueMapping.from_domain_on_another_chain(paper_apo_drm2, holo_d1.chain_id, label_seq_id_offset)  # same chain, for now, as in d1

    # create highlight script (by the spans, or just create multiple selections)
    # copy the 2 structures to 4 (paper spans vs our spans), so we can color them differently
    # select only the domains (2), and make only them visible

    sholo = superimposed_holo

    pymol_script = f"""
fetch {apo.id}
load {sholo_file_path.absolute()}

sele apo_d1, {apo.id} and chain {apo_2DA.d1.chain_id} and {get_resi_selection(apo_2DA.d1.get_spans())}
sele apo_d2, {apo.id} and chain {apo_2DA.d2.chain_id} and {get_resi_selection(apo_2DA.d2.get_spans())}
sele apo_2DA, apo_d1 or apo_d2

sele holo_d1, {sholo.id} and chain {holo_2DA.d1.chain_id} and {get_resi_selection(holo_2DA.d1.get_spans())}
sele holo_d2, {sholo.id} and chain {holo_2DA.d2.chain_id} and {get_resi_selection(holo_2DA.d2.get_spans())}
sele holo_2DA, holo_d1 or holo_d2

# copy objects, so we can color them differently
copy paper_{apo.id}, {apo.id}
copy paper_{sholo.id}, {sholo.id}

sele paper_apo_d1, paper_{apo.id} and chain {apo_2DA.d1.chain_id} and {get_resi_selection(paper_apo_drm1.get_spans())}
sele paper_apo_d2, paper_{apo.id} and chain {apo_2DA.d2.chain_id} and {get_resi_selection(paper_apo_drm2.get_spans())}
sele paper_apo_2DA, paper_apo_d1 or paper_apo_d2

sele paper_holo_d1, paper_{sholo.id} and chain {holo_2DA.d1.chain_id} and {get_resi_selection(paper_holo_drm1.get_spans())}
sele paper_holo_d2, paper_{sholo.id} and chain {holo_2DA.d2.chain_id} and {get_resi_selection(paper_holo_drm2.get_spans())}
sele paper_holo_2DA, paper_holo_d1 or paper_holo_d2

color red, apo_d1
color yellow, apo_d2
color green, holo_d1
color blue, holo_d2

color salmon, paper_apo_d1
color paleyellow, paper_apo_d2
color palegreen, paper_holo_d1
color lightblue, paper_holo_d2

# example usage: 
hide; show surface, apo_2DA
hide; show surface, paper_apo_2DA
hide; show surface, holo_2DA
hide; show surface, paper_holo_2DA

hide; show surface, apo_2DA or holo_2DA or paper_apo_2DA or paper_holo_2DA
    """

    print(pymol_script)


def main(ah_two_domain_arrangements, analyzed_domains, paper_spans):
    """

    :param ah_two_domain_arrangements:
    :param analyzed_domains: == cropped to the chains' LCS and contain only observed residues in both chains
    :return:
    """
    analyzed_domains = analyzed_domains[['full_id', 'spans', 'pdb_code', 'chain_id', 'domain_id']]

    # convert to DataFrame, so we can easily merge the arrangements with domain_info
    ah_2DAs = pd.DataFrame({'apo_2DA':  [ah[0] for ah in ah_two_domain_arrangements],
                            'holo_2DA': [ah[1] for ah in ah_two_domain_arrangements]})

    ah_2DAs = unfold_tuple_to_columns(ah_2DAs, ['apo_d1_full_id', 'apo_d2_full_id'], 'apo_2DA')
    ah_2DAs = unfold_tuple_to_columns(ah_2DAs, ['holo_d1_full_id', 'holo_d2_full_id'], 'holo_2DA')

    df = ah_2DAs.merge(analyzed_domains, left_on='apo_d1_full_id', right_on='full_id')
    df = df.merge(analyzed_domains, left_on='apo_d2_full_id', right_on='full_id', suffixes=('_apo_d1', '_apo_d2'))
    df = df.merge(analyzed_domains, left_on='holo_d1_full_id', right_on='full_id')
    df = df.merge(analyzed_domains, left_on='holo_d2_full_id', right_on='full_id', suffixes=('_holo_d1', '_holo_d2'))

    for row, _2DA_paper_spans in zip(df.itertuples(), paper_spans):
        logger.info(f'processing {row.pdb_code_apo_d1}')
        row = row._asdict()  # convert to dict, so we can easily index columns with strings

        def get_2DA(row, apo_or_holo: str):
            assert row[f'pdb_code_{apo_or_holo}_d1'] == row[f'pdb_code_{apo_or_holo}_d2'], \
                f'Domains from {apo_or_holo} 2DA should be from the same PDB structure.'

            return TwoDomainArrangement(
                row[f'pdb_code_{apo_or_holo}_d1'],
                *(DomainResidueMapping(
                    row[f'domain_id_{apo_or_holo}_{d1_or_d2}'],
                    row[f'chain_id_{apo_or_holo}_{d1_or_d2}'],
                    segment_beginnings=np.array(row[f'spans_{apo_or_holo}_{d1_or_d2}'])[:, 0].tolist(),
                    segment_ends=      np.array(row[f'spans_{apo_or_holo}_{d1_or_d2}'])[:, 1].tolist(),
                ) for d1_or_d2 in ('d1', 'd2'))
            )

        visualize_2DA(get_2DA(row, 'apo'), get_2DA(row, 'holo'), _2DA_paper_spans)


if __name__ == '__main__':
    project_logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    logging.basicConfig()

    ah_two_domain_arrangements = [
        [(('1vr6', 'A', '1vr6A01'), ('1vr6', 'A', '1vr6A02')), (('1rzm', 'A', '1vr6A01'), ('1rzm', 'A', '1vr6A02'))],
        # [(('1usg', 'A', '1usgA01'), ('1usg', 'A', '1usgA02')), (('1usi', 'A', '1usgA01'), ('1usi', 'A', '1usgA02'))],
        # [(('1za1', 'A', '1za1A01'), ('1za1', 'A', '1za1A02')), (('1q95', 'A', '1za1A01'), ('1q95', 'A', '1za1A02'))],
        # [(('1gud', 'A', '1gudA01'), ('1gud', 'A', '1gudA02')), (('1rpj', 'A', '1gudA01'), ('1rpj', 'A', '1gudA02'))],
        # [(('1jej', 'A', '1jejA01'), ('1jej', 'A', '1jejA02')), (('1jg6', 'A', '1jejA01'), ('1jg6', 'A', '1jejA02'))],
        # [(('1zol', 'A', '1zolA01'), ('1zol', 'A', '1zolA02')),	(('1o03', 'A', '1zolA01'), ('1o03', 'A', '1zolA02'))],
        # [(('1hoo', 'B', '1hooB01'), ('1hoo', 'B', '1hooB03')), (('1cg0', 'A', '1hooB01'), ('1cg0', 'A', '1hooB03'))],
        # ((('1i7n'   , 'A', '1i7nA01'), ('1i7n', 'A', '1i7nA03')), (('1i7l', 'B', '1i7nA01'), ('1i7l', 'B', '1i7nA03'))),
        # ((('5lyo', 'B', '5lyoB01'), ('5lyo', 'B', '5lyoB02')), (('4is5', 'A', '5lyoB01'), ('4is5', 'A', '5lyoB02'))),
    ]

    # probably the spans are for the apo structure (numbering is not always consistent, and is not to the LCS  - can
    # start with 4 for example)
    # not entirely sure, sometimes the numbering is from holo, see differ.txt
    paper_apo_spans_string = [
        '''(1–64)
# (65–338)''',
#         '''(1–124, 247–333)
# (125–246, 334–345)''',
#         '''(1–136, 292–310)
# (137–291)''',
#         '''(1–112, 247–288)
# (113–246)''',
#         '''(1–172, 335–351)
# (173–334)''',
#         '''(1–16, 84–221)
# (17–83)'''
#         '''(1–16, 84–221)
# (17–83)'''
    ]
    paper_spans = []
    for _2DA_string in paper_apo_spans_string:
        domains_string = re.findall(r'\([^)]+\)', _2DA_string)
        assert len(domains_string) == 2
        _2DA = [[], []]
        for d, d_string in zip(_2DA, domains_string):
            for match in re.finditer(r'([0-9]+).([0-9]+)', d_string):
                span = (int(match.group(1)), int(match.group(2)))
                d.append(span)
        paper_spans.append(_2DA)

    logger.debug(paper_spans)

    # todo aby fungovalo i s get_domains_db
    # mozna udelat novej soubor, kterej nebude fungovat s paper spans..., to bude mnohem rychlejsi, navic tohle se jeste hodit může
    domain_info_filename = Path(OUTPUT_DIR) / 'output_domains_info2021-11-13T14:21:51.339065.json'
    domain_info = pd.read_json(domain_info_filename)
    domain_info = domain_info.applymap(tuple_it)

    main(ah_two_domain_arrangements, domain_info[domain_info['type'] == 'analyzed_domain'], paper_spans)
