""" Compute LCS for all potential apo-holo pairs """
import dataclasses
import itertools
import json
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime
from difflib import SequenceMatcher
from functools import partial
from pathlib import Path
from typing import Tuple, List

import more_itertools
import numpy as np
import pandas as pd
from Bio import pairwise2

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.json_serialize import CustomJSONEncoder
from apo_holo_structure_stats.pipeline.utils.json import read_jsons_with_seqs_in, maybe_print
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_short_tasks, FutureLike, process_execute

logger = logging.getLogger(__name__)


# todo skip heterogeneity here? I skip it already in filter structures, but I wouln't know it here (as I only get the
#    get the sequence from the json, would have to add has_microheterogeneity flag or list them)
# todo mozna rovnou gzipnout json
# todo proc json.dump zere tolik pameti?
#    zere navic asi 1 - 1,5 gb, json ma pritom jen 700 MB, par vln tam ale bylo
# misto pandas jsem mohl pouzit itertools.groupby, vsechno ale musis predem .sortovat. To bude ale rychly..

"""More-
over, to avoid unresolved fragments of a structure and
misnumbering of residues, the amino acid sequences
used in this study were determined for the fragments of
structures, for which at least backbone coordinates were
available.
Subsequently, ligand-bound and ligand-free proteins
were paired at the level of 100% sequence identity"""

"""Přitom mám 30 případů, kdy substring není maximalní. To může být způsobeno:
- brali jen resolved fragments
- porovnávám 3letter cody, nevim, co dělali oni

na jednom přikladu to tak je, ze zacatek berou jen observed """


@dataclasses.dataclass
class LCSResult:
    # lcs: List[str]  # three letter codes
    length: int  # substring length
    i1: int
    i2: int
    # i1_seq_num: int  # first sequence (apo) start label_seq_id
    # i2_seq_num: int  # second sequence (holo) start label_seq_id
    mismatches: int  # leading_mismatches + trailing_mismatches
    leading_mismatches: int  # mismatches before LCS
    trailing_mismatches: int  # mismatches after LCS


# todo multiprocessing logger -- need to prepend the current processed chain...
# or logger for each process and prepend process id (that would work

MIN_SUBSTRING_LENGTH_RATIO = 0.9


def get_longest_common_polypeptide(
        apo_seq: List[str], holo_seq: List[str],  # in 3-letter codes
):
    """ Substring """

    # find longest common substring, works even with the 3-letter codes (sequence of strings), e.g. MET and FME will get treated differently
    l1, l2 = len(apo_seq), len(holo_seq)
    i1, i2, length = SequenceMatcher(a=apo_seq, b=holo_seq, autojunk=False).find_longest_match(0, l1, 0, l2)

    #todo tohle kontrolovat neni tolik informativni, me zajima, jestli tam je mismatch, ale muze se stat, ze jsou jenom
    #  posunuty (s1 ma leading cast a s2 trailing),
    # chci vyprintovat pocet mismatchu (zacatek a konec),
    # jednoduše zjistim z i1, i2 (min z nich pocet mismatchu na zacatku)
    # min (l1 - i1+length, l2 - i2 + length) pocet mismatchu na konci
    substring_length_ratio = length / min(l1, l2)

    leading_mismatches = min(i1, i2)
    trailing_mismatches = min(l1 - (i1 + length), l2 - (i2 + length))
    mismatches = leading_mismatches + trailing_mismatches

    # todo obcas jsou jenom dva mismatche v aligmentu (ne substring), ma cenu prozkoumavat to?
    #
    # if mismatches > 0:
    #     # alignment = next(aligner.align(apo_residue_codes, holo_residue_codes))
    #     alignment = pairwise2.align.globalms(apo_seq, holo_seq,
    #                                          1, 0,  # match, mismatch
    #                                          -.5, -.1, # gap open, ext
    #                                          penalize_end_gaps=False, gap_char=['-'], one_alignment_only=True)[0]
    #     logger.debug('Sequences differ, alignment:')
    #     # logger.info(f'\n{alignment}')
    #     logger.debug(f'\n{pairwise2.format_alignment(*alignment)}')
    #     logger.info(f'found {mismatches} mismatches ({leading_mismatches} substring leading, '
    #                     f'{trailing_mismatches} trailing)')
    #     logger.info(f'substring length: {substring_length_ratio:.3f} < {MIN_SUBSTRING_LENGTH_RATIO}')

    # # crop polypeptides to longest common substring
    # apo_common_seq = dict(itertools.islice(apo_seq.items(), i1, i1+length))
    # holo_common_seq = dict(itertools.islice(holo_seq.items(), i2, i2+length))
    #
    # # return the sequences (values will be same, but not the keys, label_seq_ids, they might be offset, or depending on its definition
    # # (which I find ambiguous), in a more complicated relation)

    return LCSResult(
        # apo_seq[i1: i1 + length],  # would take unnecessary memory (one could easily get the sequence with length and i)
        length,
        i1,
        i2,
        # i1 + next(iter(apo_seq.keys())),  # todo tohle nemám, pokud jenom ukládám sekvenci a ne start label_seq_id
        # i2 + next(iter(holo_seq.keys())),
        mismatches,
        leading_mismatches,
        trailing_mismatches,
    )


# cannot be in a lambda, as I use multiprocessing (not picklable)
def fn_wrapper_unpack_args(fn, args):
    return fn(*args)


def compute_lcs(apo_seq: List[str], holo_seq: List[str],
                apo_chain: Tuple[str, str], holo_chain: Tuple[str, str], group_label):
    logger.info(f'processing {apo_chain} {holo_chain}, from {group_label}')

    # run lcs
    lcs_result = get_longest_common_polypeptide(apo_seq, holo_seq)
    return lcs_result

#
# def get_pairs_in_group(structures_metadata_group, group_label=None):
#     df = structures_metadata_group
#     apo = df[~df.is_holo]
#     holo = df[df.is_holo]
#
#     logger.info(f'{group_label}: apo: {len(apo)}, holo: {len(holo)}')
#
#     # product
#     potential_pairs = apo.merge(holo, how='cross', suffixes=('_apo', '_holo'))
#
#     # return futures
#     for row in potential_pairs.itertuples():
#         apo_chain = (row.pdb_code_apo, row.chain_id_apo)
#         holo_chain = (row.pdb_code_holo, row.chain_id_holo)
#
#         apo_seq = row.sequence_apo
#         holo_seq = row.sequence_holo
#
#         # yield a tuple of metadata and args for `get_longest_common_polypeptide`
#         yield (apo_chain, holo_chain), (apo_seq, holo_seq)
#

# def get_pairs_in_group(structures_metadata, group_indices, group_label=None):
#     df = structures_metadata.iloc[group_indices]
#     group_is_holo = structures_metadata.is_holo.iloc[group_indices]
#     apo_indices = group_indices[~group_is_holo]
#     holo_indices = group_indices[group_is_holo]
#
#     logger.info(f'{group_label}: apo: {len(apo_indices)}, holo: {len(holo_indices)}')
#
#     # product apo*holo
#     for apo_index in apo_indices:
#         apo_row = structures_metadata.iloc[apo_index]
#         apo_chain = (apo_row.pdb_code, apo_row.chain_id)
#
#         for holo_index in holo_indices:
#             holo_row = structures_metadata.iloc[holo_index]
#             holo_chain = (holo_row.pdb_code, holo_row.chain_id)
#
#             # yield a tuple of metadata and args for `get_longest_common_polypeptide`
#             yield (apo_chain, holo_chain), (apo_row.sequence, holo_row.sequence)

def get_pairs_in_group(structures_metadata, group_indices, group_label=None):
    df = structures_metadata.iloc[group_indices]
    apo = df[~df.is_holo]
    holo = df[df.is_holo]

    logger.info(f'{group_label}: apo: {len(apo)}, holo: {len(holo)}')

    holo_tuples = list(holo.itertuples())  # put in a list, so we call itertuples() only once (not in the loop)

    # product apo*holo
    for apo_row in apo.itertuples():
        apo_chain = (apo_row.pdb_code, apo_row.chain_id)

        for holo_row in holo_tuples:
            holo_chain = (holo_row.pdb_code, holo_row.chain_id)

            # yield a tuple of metadata and args for `get_longest_common_polypeptide`
            yield (apo_chain, holo_chain), (apo_row.sequence, holo_row.sequence)


def make_pairs_with_lcs(structures_metadata, workers):    # groupby isoform
    # product apo holo
    # do following in multiple processes:
    #   find LCS
    # write LCS to output (if no mismatches?), or all?

    groups = structures_metadata.groupby('uniprotkb_id')  # or ísoform (maybe specify in args)
    # structures as files or codes, so that should be handled in `chains_for_uniprot_ids`? Or somehow with a join in filter_structures?
    # that could be done..
    def get_pairs():
        for uniprot_id, group_indices in groups.indices.items():
            for pair in get_pairs_in_group(structures_metadata, group_indices, uniprot_id):
                yield uniprot_id, pair

    uniprot_ids, pairs = more_itertools.unzip(get_pairs())
    pair_ids, lcs__args = more_itertools.unzip(pairs)

    i = 0
    print(datetime.now())

    # for uniprot_id, (apo_chain, holo_chain), lcs_future in zip(uniprot_ids, pair_ids, lcs_futures):
    for uniprot_id, (apo_chain, holo_chain), args in zip(uniprot_ids, pair_ids, lcs__args):
        lcs_future = FutureLike(process_execute(get_longest_common_polypeptide, *args))
        i += 1
        if i % 100 == 0:
            # if i % 100 == 0:
            maybe_print(False, f'\r{i}', end='')

        try:
            logger.info(f'getting result of {apo_chain} {holo_chain}, from {uniprot_id}')

            yield {
                'pdb_code_apo': apo_chain[0],  # todo could rename the variables to more general chain1 /c1, c2...
                'chain_id_apo': apo_chain[1],
                'pdb_code_holo': holo_chain[0],
                'chain_id_holo': holo_chain[1],
                'lcs_result': lcs_future.result(),
            }

        except Exception as e:
            logger.exception('compute_lcs failed with: ')

    print(datetime.now())


def make_pairs_with_lcs_old(structures_metadata, workers):
    # groupby isoform
    # product apo holo
    # do following in multiple processes:
    # find LCS
    # write LCS to output (if no mismatches?), or all?

    groups = structures_metadata.groupby('uniprotkb_id')  # or ísoform (maybe specify in args)
    # structures as files or codes, so that should be handled in `chains_for_uniprot_ids`? Or somehow with a join in filter_structures?
    # that could be done..
    print(list(groups.indices.items())[:10])
    quit()

    # with ProcessPoolExecutor(workers) as executor:
    with ThreadPoolExecutor(1) as executor:
        def get_pairs():
            for uniprot_id, group in groups:
                for pair in get_pairs_in_group(group, uniprot_id):
                    yield uniprot_id, pair
                # yield from get_pairs_in_group(group, uniprot_id)

        uniprot_ids, pairs = more_itertools.unzip(get_pairs())
        pair_ids, lcs__args = more_itertools.unzip(pairs)

        # lcs_futures = submit_short_tasks(executor, 40 * workers, 100,
        # lcs_futures = submit_short_tasks(executor, 4 * workers, 10,
        #                                  partial(fn_wrapper_unpack_args, get_longest_common_polypeptide),
        #                                  lcs__args)

        # PYPY, jen iter
        # trva dlouho jenom proiterovat, obcas jsou tam takovy zachvevy - to je asi kvuli groupby, ze je vytvari df
        # pro kazdou skupinu? Jasne, to cross product zejo atd.. Mohl bych to nechat v pythonu ciste, avt to slo rychleji
        # jenom iterovani trva: 15K párů za 4 min (ale u tech velkejch, kt. je nejvic je to hodne rychly)
        # nakonec to dalo ty 3M za 14 min, 5,5 GB mem (vysledky asi ), 7 GB kdyz zapisuje na disk, hodne pomaly
        # json ma nakonec "jenom" 670 MB??
        i = 0
        print(datetime.now())
        # for uniprot_id, (apo_chain, holo_chain), lcs_args in zip(uniprot_ids, pair_ids, lcs__args):
        #     i += 1
        #     if i % 10 == 0:
        #         # if i % 100 == 0:
        #         maybe_print(False, f'\r{i}', end='')
        #
        #     yield {
        #         'pdb_code_apo': apo_chain[0],
        #         'chain_id_apo': apo_chain[1],
        #         'pdb_code_holo': holo_chain[0],
        #         'chain_id_holo': holo_chain[1],
        #         'lcs_result': LCSResult(i, i//2, i//2+1, 0,0,0),
        #     }
        # return
        # fn = partial(fn_wrapper_unpack_args, get_longest_common_polypeptide)
        # lcs_futures = (FutureLike(process_execute(fn, args)) for args in lcs__args)

        # for uniprot_id, (apo_chain, holo_chain), lcs_future in zip(uniprot_ids, pair_ids, lcs_futures):
        for uniprot_id, (apo_chain, holo_chain), args in zip(uniprot_ids, pair_ids, lcs__args):
            lcs_future = FutureLike(process_execute(get_longest_common_polypeptide, *args))
            i += 1
            if i % 10 == 0:
                # if i % 100 == 0:
                maybe_print(False, f'\r{i}', end='')

            try:
                logger.info(f'getting result of {apo_chain} {holo_chain}, from {uniprot_id}')


                yield {
                    'pdb_code_apo': apo_chain[0],
                    'chain_id_apo': apo_chain[1],
                    'pdb_code_holo': holo_chain[0],
                    'chain_id_holo': holo_chain[1],
                    'lcs_result': lcs_future.result(),
                }

            except Exception as e:
                logger.exception('compute_lcs failed with: ')



def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=1, help='number of threads for concurrent API requests')
    parser.add_argument('structures_json', type=Path, help='File needs to contain list of objects with pdb_code, chain_id, isoform '
                                                'and is_holo flag')
    parser.add_argument('output_file', type=Path, help='writes apo-holo pairs in json')
    add_loglevel_args(parser)
    args = parser.parse_args()
    # todo combine this and put in logs (I only use this in scripts anyway)
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    structures_metadata = read_jsons_with_seqs_in(str(args.structures_json), quiet=False)#args.loglevel > logging.INFO)
    pairs = list(make_pairs_with_lcs(structures_metadata, args.workers))

    with args.output_file.open('w') as f:
        json.dump(pairs, f, cls=CustomJSONEncoder)


if __name__ == '__main__':
    main()


def load_pairs_json(pairs_json: str):
    def convert_lcs_result(items):
        try:
            items['lcs_result'] = LCSResult(**items['lcs_result'])
        except KeyError:
            pass

        return items

    with open(pairs_json) as f:
        return pd.DataFrame.from_records(json.load(f, object_hook=convert_lcs_result))


def pairs_without_mismatches(potential_pairs: pd.DataFrame):
    return potential_pairs[0 == np.array([lcs_result.mismatches for lcs_result in potential_pairs.lcs_result])]
