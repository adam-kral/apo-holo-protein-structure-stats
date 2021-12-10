""" Compute LCS for all potential apo-holo pairs """
import concurrent
import dataclasses
import itertools
import json
import logging
from concurrent.futures import ProcessPoolExecutor
from difflib import SequenceMatcher
from pathlib import Path
from typing import Tuple, List, Dict

import pandas as pd
from Bio import pairwise2

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.json_serialize import CustomJSONEncoder
from apo_holo_structure_stats.input.download import parse_mmcif
from apo_holo_structure_stats.pipeline.log import add_loglevel_args

logger = logging.getLogger(__name__)


# todo skip heterogeneity here
# todo skip mismatches here
# todo save LCS as label_seq_id begin in s1 and s2, plus the list? (or joined string) of three letter codes

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

# todo combine filter_structures and is_holo to one script? (But still would want to know when something failed,
#  is holo shouldn't fail though, it has no filtering, only classification)


@dataclasses.dataclass
class LCSResult:
    lcs: List[str]  # three letter codes
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
        apo_seq: Dict[int, str], holo_seq: Dict[int, str],  # in 3-letter codes
):
    """ Substring """

    # find longest common substring, works even with the 3-letter codes (sequence of strings), e.g. MET and FME will get treated differently
    apo_residue_codes = list(apo_seq.values())
    holo_residue_codes = list(holo_seq.values())
    l1, l2 = len(apo_seq), len(holo_seq)
    i1, i2, length = SequenceMatcher(a=apo_residue_codes, b=holo_residue_codes, autojunk=False) \
        .find_longest_match(0, l1, 0, l2)

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

    if mismatches > 0:
        # alignment = next(aligner.align(apo_residue_codes, holo_residue_codes))
        alignment = pairwise2.align.globalms(apo_residue_codes,holo_residue_codes,
                                             1, 0,  # match, mismatch
                                             -.5, -.1, # gap open, ext
                                             penalize_end_gaps=False, gap_char=['-'], one_alignment_only=True)[0]
        logger.debug('Sequences differ, alignment:')
        # logger.info(f'\n{alignment}')
        logger.debug(f'\n{pairwise2.format_alignment(*alignment)}')
        logger.warning(f'found {mismatches} mismatches ({leading_mismatches} substring leading, '
                        f'{trailing_mismatches} trailing)')
        logger.warning(f'substring length: {substring_length_ratio:.3f} < {MIN_SUBSTRING_LENGTH_RATIO}')

    # # crop polypeptides to longest common substring
    # apo_common_seq = dict(itertools.islice(apo_seq.items(), i1, i1+length))
    # holo_common_seq = dict(itertools.islice(holo_seq.items(), i2, i2+length))
    #
    # # return the sequences (values will be same, but not the keys, label_seq_ids, they might be offset, or depending on its definition
    # # (which I find ambiguous), in a more complicated relation)

    return LCSResult(
        apo_residue_codes[i1: i1 + length],
        length,
        i1,
        i2,
        # i1 + next(iter(apo_seq.keys())),  # todo tohle nemám, pokud jenom ukládám sekvenci a ne start label_seq_id
        # i2 + next(iter(holo_seq.keys())),
        mismatches,
        leading_mismatches,
        trailing_mismatches,
    )


def compute_lcs(apo_chain: Tuple[str, str], holo_chain: Tuple[str, str], group_label=None):
    apo_pdb_code, apo_chain_id = apo_chain
    holo_pdb_code, holo_chain_id = holo_chain

    logger.info(f'processing {apo_chain} {holo_chain}, from {group_label}')

    # load the structure from file
    apo, (apo_residue_id_mappings, apo_poly_seqs) = parse_mmcif(apo_pdb_code)
    holo, (holo_residue_id_mappings, holo_poly_seqs) = parse_mmcif(holo_pdb_code)

    # get the first model (s[0]), (in x-ray structures there is probably a single model)
    apo, apo_residue_id_mappings = map(lambda s: s[0], (apo, apo_residue_id_mappings))
    holo, holo_residue_id_mappings = map(lambda s: s[0], (holo, holo_residue_id_mappings))

    # todo tady spadne microheterogeneity, asi lepsi, kdyby ve filter structures?
    apo_mapping = apo_residue_id_mappings[apo_chain_id]
    holo_mapping = holo_residue_id_mappings[holo_chain_id]

    # run lcs
    lcs_result = get_longest_common_polypeptide(apo_poly_seqs[apo_mapping.entity_poly_id],
                                                holo_poly_seqs[holo_mapping.entity_poly_id])
    return lcs_result


def make_pairs_with_lcs_in_group(structures_metadata_group, executor, group_label=None):
    df = structures_metadata_group
    apo = df[~df.is_holo]
    holo = df[df.is_holo]

    # product
    potential_pairs = apo.merge(holo, how='cross', suffixes=('_apo', '_holo'))

    # return futures
    for row in potential_pairs.itertuples():
        apo_chain = (row.pdb_code_apo, row.chain_id_apo)
        holo_chain = (row.pdb_code_holo, row.chain_id_holo)

        future = executor.submit(compute_lcs, apo_chain, holo_chain, group_label)

        yield apo_chain, holo_chain, future


def make_pairs_with_lcs(structures_metadata, workers):
    # groupby isoform
        # product apo holo
            # do following in multiple processes:
            # find LCS
            # write LCS to output (if no mismatches?), or all?

    groups = structures_metadata.groupby('isoform')

    chains_and_futures = []

    with ProcessPoolExecutor(workers) as executor:
        for uniprot_id, group in groups:
            chains_and_futures.extend(make_pairs_with_lcs_in_group(group, executor, uniprot_id))

    for apo_chain, holo_chain, lcs_future in chains_and_futures:
        try:
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

    structures_metadata = pd.read_json(args.structures_json).iloc[:100]
    pairs = list(make_pairs_with_lcs(structures_metadata, args.workers))

    with args.output_file.open('w') as f:
        json.dump(pairs, f, cls=CustomJSONEncoder)


if __name__ == '__main__':
    main()
