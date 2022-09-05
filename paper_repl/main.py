import concurrent.futures
import itertools
import json
import logging
from datetime import datetime
from difflib import SequenceMatcher
from multiprocessing import Manager
from pathlib import Path
from typing import List, Dict

import pandas as pd
from Bio import pairwise2

from Bio.PDB.Chain import Chain
from Bio.pairwise2 import format_alignment

from apo_holo_structure_stats.core.dataclasses import ChainResidueData, DomainResidueData, ChainResidues, \
    DomainResidueMapping, DomainResidues
from apo_holo_structure_stats.core.analysesinstances import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.input.download import APIException, parse_mmcif
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds, ResidueId
from apo_holo_structure_stats.pipeline.run_analyses import AnalysisHandler, JSONAnalysisSerializer, \
    ConcurrentJSONAnalysisSerializer, get_observed_residues

from .paper_dataset_utils import get_chain_by_chain_code

OUTPUT_DIR = 'output'


def get_longest_common_polypeptide(
        apo_seq: Dict[int, str], holo_seq: Dict[int, str],  # in 3-letter codes
    ):
    """ Substring """

    # find longest common substring, works even with the 3-letter codes (sequence of strings), e.g. MET and FME will get treated differently
    apo_residue_codes = list(apo_seq.values())
    holo_residue_codes = list(holo_seq.values())
    l1, l2 = len(apo_seq), len(holo_seq)
    i1, i2, length = SequenceMatcher(a=apo_residue_codes, b=holo_residue_codes, autojunk=False)\
        .find_longest_match(0, l1, 0, l2)

    #todo tohle kontrolovat neni tolik informativni, me zajima, jestli tam je mismatch, ale muze se stat, ze jsou jenom
    #  posunuty (s1 ma leading cast a s2 trailing),
    # chci vyprintovat pocet mismatchu (zacatek a konec),
    # jednoduše zjistim z i1, i2 (min z nich pocet mismatchu na zacatku)
    # min (l1 - i1+length, l2 - i2 + length) pocet mismatchu na konci
    substring_length_ratio = length / min(l1, l2)
    logging.info(f'substring to original ratio: {substring_length_ratio}')

    leading_mismatches = min(i1, i2)
    trailing_mismatches = min(l1 - (i1 + length), l2 - (i2 + length))
    mismatches = leading_mismatches + trailing_mismatches

    if mismatches > 0:
        # alignment = next(aligner.align(apo_residue_codes, holo_residue_codes))
        alignment = pairwise2.align.globalms(apo_residue_codes,holo_residue_codes,
                                             1, 0,  # match, mismatch
                                             -.5, -.1, # gap open, ext
                                             penalize_end_gaps=False, gap_char=['-'], one_alignment_only=True)[0]
        logging.info('Sequences differ, alignment:')
        # logging.info(f'\n{alignment}')
        logging.info(f'\n{format_alignment(*alignment)}')
        logging.warning(f'found {mismatches} mismatches ({leading_mismatches} substring leading, '
                        f'{trailing_mismatches} trailing)')
        logging.warning(f'substring length: {substring_length_ratio:.3f} < {MIN_SUBSTRING_LENGTH_RATIO}')

    # crop polypeptides to longest common substring
    apo_common_seq = dict(itertools.islice(apo_seq.items(), i1, i1+length))
    holo_common_seq = dict(itertools.islice(holo_seq.items(), i2, i2+length))

    # return the sequences (values will be same, but not the keys, label_seq_ids, they might be offset, or depending on its definition
    # (which I find ambiguous), in a more complicated relation)
    return apo_common_seq, holo_common_seq


MIN_SUBSTRING_LENGTH_RATIO = 0.90


if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)


    def run_apo_analyses():
        df = pd.read_csv('apo_holo.dat', delimiter=r'\s+', comment='#', header=None,
                         names=('apo', 'holo', 'domain_count', 'ligand_codes'), dtype={'domain_count': int})

        start_datetime = datetime.now()
        analyses_output_fpath = Path(OUTPUT_DIR) / f'output_apo_holo_{start_datetime.isoformat()}.json'
        domains_info_fpath = Path(OUTPUT_DIR) / f'output_domains_info_{start_datetime.isoformat()}.json'

        with Manager() as multiprocessing_manager:

            serializer = ConcurrentJSONAnalysisSerializer(analyses_output_fpath, multiprocessing_manager)
            domains_info = multiprocessing_manager.list()

            found = False

            # with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
                for index, row in df.iterrows():
                    # if row.apo[:4] == '1cgj':
                    #     continue  # chymotrypsinogen x chymotrypsin + obojí má ligand... (to 'apo' má 53 aa inhibitor)

                    # if row.apo[:4] == '1ikp':
                    #     found = True

                    # if row.apo[:4] not in ('2d6l', ):
                    #     continue
                    #
                    # if not found:
                    #     continue

                    future = executor.submit(process_pair,
                        s1_pdb_code=row.apo[:4],
                        s2_pdb_code=row.holo[:4],
                        s1_paper_chain_code=row.apo[4:],
                        s2_paper_chain_code=row.holo[4:],
                        serializer=serializer,
                        domains_info=domains_info,
                    )

            serializer.dump_data()
            with domains_info_fpath.open('w') as f:
                json.dump(list(domains_info), f)

            print(start_datetime.isoformat())
            print(datetime.now().isoformat())

    def run_holo_analyses():
        df = pd.read_csv('holo.dat', delimiter=r'\s+', comment='#', header=None,
                         names=('apo', 'holo', 'domain_count', 'ligand_codes'), dtype={'domain_count': int})


        serializer = JSONAnalysisSerializer('output_holo_all.json')


        found = False
        for index, row in df.iterrows():
            # if row.apo[:4] == '1cgj':
            #     continue  # chymotrypsinogen x chymotrypsin + obojí má ligand... (to 'apo' má 53 aa inhibitor)
            #
            # if row.apo[:4] == '1bq8':
            #     found = True
            #
            # if not found:
            #     continue

            # if row.apo[:4] != '1bq8':
            #     continue

            logging.info(f'{row.apo[:4]}, {row.holo[:4]}')

            apo = parse_mmcif(row.apo[:4])
            holo = parse_mmcif(row.holo[:4])

            apo_chain = get_chain_by_chain_code(apo, row.apo[4:])
            holo_chain = get_chain_by_chain_code(holo, row.holo[4:])

            try:
                compare_chains(apo_chain, holo_chain,
                               [get_rmsd],
                               [compare_ss],
                               [get_rmsd],
                               [compare_ss],
                               [get_hinge_angle],
                               serializer)
            except Exception as e:
                logging.exception('compare chains failed with: ')
                # raise

        pass


    run_apo_analyses()
    # run_holo_analyses()


