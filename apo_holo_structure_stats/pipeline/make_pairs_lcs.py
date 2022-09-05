""" Pair chains for subsequent structural analyses of the pairs.

Default impl. computes the longest common substring for all potential apo-holo pairs within a uniprot accession.
"""
import argparse
import dataclasses
import glob
import itertools
import json
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, Executor
from datetime import datetime
from difflib import SequenceMatcher
from functools import partial
from pathlib import Path
from typing import Tuple, List, Iterable, NamedTuple

import more_itertools
import numpy as np
import pandas as pd

from apo_holo_structure_stats import project_logger, settings
from apo_holo_structure_stats.core.json_serialize import CustomJSONEncoder
from apo_holo_structure_stats.pipeline.utils.json import read_jsons_with_seqs, maybe_print
from apo_holo_structure_stats.pipeline.utils.log import get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_short_tasks
from apo_holo_structure_stats.settings import Settings

logger = logging.getLogger(__name__)


# todo mozna rovnou gzipnout json


@dataclasses.dataclass
class LCSResult:
    length: int  # substring length
    i1: int
    i2: int
    mismatches: int  # leading_mismatches + trailing_mismatches
    leading_mismatches: int  # mismatches before LCS
    trailing_mismatches: int  # mismatches after LCS


def get_longest_common_polypeptide(
        apo_seq: List[str], holo_seq: List[str],  # in 3-letter codes
):
    """ Find the longest common substring, works even with the 3-letter codes (sequence of strings),
    e.g. MET and FME will get treated differently"""

    l1, l2 = len(apo_seq), len(holo_seq)
    i1, i2, length = SequenceMatcher(a=apo_seq, b=holo_seq, autojunk=False).find_longest_match(0, l1, 0, l2)

    leading_mismatches = min(i1, i2)
    trailing_mismatches = min(l1 - (i1 + length), l2 - (i2 + length))
    mismatches = leading_mismatches + trailing_mismatches

    return LCSResult(
        length,
        i1,
        i2,
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

# was too slow, as mentioned below (cross), at least in pypy because of c pandas
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

# was this faster than the current impl? I don't remember. Probably? But iloc creates a Series.
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


class Matchmaker:
    """ Pairs equivalent chains for subsequent analyses. Outputs pairs with metadata for each pair.

    By default, the longest common substring is computed and mismatches preceding and trailing the lcs are counted.
    When equivalent chains would be those without the mismatches.

    By default, all chains within a uniprot accession are paired together (all combinations).

    See make_pairs for more details.
    """
    def __init__(self, workers: int):
        self.workers = workers

    def product_pairs(self, df1, df2) -> Iterable[Tuple[Tuple, Tuple]]:
        return itertools.product(df1.itertuples(), df2.itertuples())

    def combination_pairs(self, df) -> Iterable[Tuple[Tuple, Tuple]]:
        return itertools.combinations(df.itertuples(), 2)

    def get_pairs_in_group(self, chain_group, group_label=None) -> Iterable[Tuple[Tuple, Tuple]]:
        return self.combination_pairs(chain_group)

    def get_potential_pairs(self, chain_metadata):
        """
        :param chain_metadata: pandas.DataFrame
        :return: generator of tuples of uniprot_id, (c1_namedtuple, c2_namedtuple)
        """
        groups = chain_metadata.groupby('uniprotkb_id')
        for uniprot_id, group_indices in groups.indices.items():
            for pair in self.get_pairs_in_group(chain_metadata.iloc[group_indices], uniprot_id):
                yield uniprot_id, pair

    def get_pair_dict(self, uniprot_id: str, c1: NamedTuple, c2, lcs_result: LCSResult):
        return {
            'pdb_code_c1': c1.pdb_code,
            'chain_id_c1': c1.chain_id,
            'pdb_code_c2': c2.pdb_code,
            'chain_id_c2': c2.chain_id,
            'lcs_result': lcs_result,
        }

    def make_pairs(self, chain_metadata: pd.DataFrame) -> Iterable[dict]:
        """ Pairs equivalent chains for subsequent analyses. Returns pairs with metadata for each pair.

        Longest common substring is computed and mismatches preceding and trailing the lcs are counted.
        When equivalent chains would be those without the mismatches.

        To reduce the number of pairs the lcs is computed for, we group the chains by uniprotkb_id and pair them
        only within the same group.

        :param chain_metadata: pandas.DataFrame
            - columns: pdb_code, chain_id, uniprotkb_id, sequence
                - pdb_code: string
                - chain_id: string
                - uniprotkb_id: string, primary uniprot accession
                - sequence: list of 3-letter amino acid codes

        :return: generator of dictionaries of form as in self.get_pair_dict.
        """
        try:
            uniprot_ids, pairs = more_itertools.unzip(self.get_potential_pairs(chain_metadata))
        except ValueError:
            # if iterable get_pairs() is empty, unpacking of more_itertools.unzip() raises this exception
            # no time to fix that in a reasonable way
            return

        pairs, pairs_copy = itertools.tee(pairs)  # pairs are tuples of apo_tuple, holo_tuple

        # parallelize the computation of LCS
        with ProcessPoolExecutor(self.workers) as executor:
            lcs_futures = submit_short_tasks(executor, 4 * self.workers, 100,
                                             partial(fn_wrapper_unpack_args, get_longest_common_polypeptide),
                                             ((c1.sequence, c2.sequence) for c1, c2 in pairs))

            i = 0
            for uniprot_id, (c1, c2), lcs_future in zip(uniprot_ids, pairs_copy, lcs_futures):
                i += 1
                if i % 100 == 0:
                    maybe_print(logger.level <= logging.INFO, f'\r{i}', end='')
                    logger.info(f'the last hundredth processed pair: {(c1.pdb_code, c1.chain_id)} '
                                f'{(c2.pdb_code, c2.chain_id)}, from {uniprot_id}')

                try:
                    yield self.get_pair_dict(uniprot_id, c1, c2, lcs_future.result())
                except Exception as e:
                    logger.exception('compute_lcs failed with: ')


class ApoHoloMatchmaker(Matchmaker):
    """ Makes pairs consisting of an apo and holo chain, given the chain metadata.

    In addition to the requirements in Matchmaker.make_pairs, this class also requires field `is_holo` to be
    present in the chain metadata. """
    def get_pairs_in_group(self, chain_group, group_label=None):
        df = chain_group
        apo = df[~df.is_holo]
        holo = df[df.is_holo]
        logger.info(f'{group_label}: apo: {len(apo)}, holo: {len(holo)}')

        return self.product_pairs(apo, holo)

    def get_pair_dict(self, uniprot_id, c1, c2, lcs_result):
        return {
            'pdb_code_apo': c1.pdb_code,
            'chain_id_apo': c1.chain_id,
            'pdb_code_holo': c2.pdb_code,
            'chain_id_holo': c2.chain_id,
            'lcs_result': lcs_result,
        }


# PYPY, jen iter
# trva dlouho jenom proiterovat argumenty, obcas jsou tam takovy zachvevy - to je asi kvuli groupby, ze je vytvari df
# pro kazdou skupinu? Jasne, to cross product zejo atd.. Mohl bych to nechat v pythonu ciste, avt to slo rychleji
# jenom iterovani trva: 15K párů za 4 min (ale u tech velkejch, kt. je nejvic je to hodne rychly)
# nakonec to dalo ty 3M za 14 min, 5,5 GB mem (vysledky asi ), 7 GB kdyz zapisuje na disk, hodne pomaly
# json ma nakonec "jenom" 670 MB??


parser = get_argument_parser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=__doc__,
                             epilog='''\
To modify the behavior of this script, set Settings.MAKE_PAIRS_CLASS to your subclass of Matchmaker.

Creates JSON with records for each potential pair (within a uniprot accesion) with fields:
- pdb_code_apo, chain_id_apo, pdb_code_holo, chain_id_holo, lcs_result (see LCSResult class)
- use `load_pairs_json` to load the JSON into a pandas.DataFrame and `pairs_without_mismatches` to filter out
potential pairs with mismatches leading or trailing the LCS.

''')
parser.add_argument('--workers', type=int, default=1, help='Number of workers to use when computing LCS.')
parser.add_argument('structures_json', type=Path, help='JSON containing of records with pdb_code, '
                                                       'chain_id, uniprotkb_id, and is_holo flag')
parser.add_argument('output_file', type=Path, help='writes apo-holo pairs in json')

def main():
    args = parser.parse_args()
    # todo combine this and put in logs (I only use this in scripts anyway)
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()

    json_files = glob.glob(str(args.structures_json))
    chain_metadata = read_jsons_with_seqs(json_files, quiet=False)#args.loglevel > logging.INFO)

    matchmaker_class = settings.load_class_from_str(Settings.MAKE_PAIRS_CLASS)
    matchmaker: Matchmaker = matchmaker_class(args.workers)
    pairs = list(matchmaker.make_pairs(chain_metadata))

    logger.info(f'Computed LCS for {len(pairs)} pairs.')
    logger.info(f'Writing to {args.output_file}...')

    with args.output_file.open('w') as f:
        json.dump(pairs, f, cls=CustomJSONEncoder)

    logger.info(f'Done.')

""" """
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
