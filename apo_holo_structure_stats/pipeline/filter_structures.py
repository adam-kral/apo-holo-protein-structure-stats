#!/usr/bin/env python3
import concurrent
import itertools
import json
import os
import warnings
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path
from typing import Set

import pandas as pd
from Bio.PDB import PDBList

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analyses import GetChains, IsHolo
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds
from apo_holo_structure_stats.input.download import find_or_download_structure, parse_mmcif
from apo_holo_structure_stats.pipeline.download_structures import download_structures
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args, get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks
from apo_holo_structure_stats.settings import Settings
from apo_holo_structure_stats.core.dataclasses import ChainResidues

logger = logging.getLogger(__name__)


def get_resolution(mmcif_dict):
    """ Basically the same as in BioPython, but fixes it for cryo-em structures (on which BioPython returns None, because
    they often or always have `.` in _refine.ls_d_res_high.)

    https://github.com/biopython/biopython/pull/3229#issuecomment-678608377
    """
    res_keys = [
        "_refine.ls_d_res_high",
        "_refine_hist.d_res_high",
        "_em_3d_reconstruction.resolution",
    ]

    resolution = None

    for key in res_keys:
        try:
            val = mmcif_dict[key]
            item = val[0]
        except (KeyError, IndexError):
            continue
        if item not in ('?', '.'):
            resolution = item

    if resolution is not None:
        try:
            resolution = float(resolution)
        except ValueError:
            pass

    return resolution


def structure_meets_our_criteria(s, s_header, mmcif_dict, get_chains: GetChains):
    # todo neni s structure, ale model!
    """ decides if structure meets criteria for resolution, etc. """

    # resolution = s_header['resolution']  # this does not work with elmi
    resolution = get_resolution(mmcif_dict)
    s_id = s.get_parent().id

    # skip low resolution
    if resolution and resolution > Settings.MIN_STRUCTURE_RESOLUTION:
        logger.info(f'skipping structure {s_id}: resolution ({resolution}) does not meet the limit of '
                    f'{Settings.MIN_STRUCTURE_RESOLUTION}')
        return False
    #
    # # skip non-xray. Done in the original paper. Try without it. Or then, specify in output dataset which experimental
    # #   method (or one could integrate the info themselves)
    # if mmcif_dict['_exptl.method'][0] != 'X-RAY DIFFRACTION':
    #     logger.info(f'skipping structure {s_id}: exp. method `{mmcif_dict["_exptl.method"]}` is not X-RAY DIFFRACTION')
    #     return False

    # this is bad, this also should be in the output files, so we can easily know how many structures/chains were
    # discarded and why?
    # we already parse them, at least the mmcif
    # probably I didn't have it in the file from the beginning, as this is an attribute of the structure, not the chains...
    # so what? I don't need to create a flat file (or I could save it in long form). But then, it would be hard to read
    # it with pandas (json normalize, not exactly intuitive use -  have to list all outer columns manually)
    # no time for that now - ignore
    if not resolution:
        # for those methods, the resolution should be obtained:
        # assert mmcif_dict["_exptl.method"][0] not in ('X-RAY DIFFRACTION',  'ELECTRON MICROSCOPY')
        # skip NMR structures and others
        logger.info(f'skipping structure {s_id}: exp. method `{mmcif_dict["_exptl.method"]}`, could not determine the resolution')
        return False

    # skip DNA/RNA complexes
    try:
        if any(poly_type in mmcif_dict['_entity_poly.type'] for poly_type in ('polydeoxyribonucleotide', 'polyribonucleotide')):
            logger.info(f'skipping structure {s_id}: no interest in structures with RNA or DNA chains')
            return False
    except KeyError:
        # _entity_poly.type:: Required in PDB entries no; Used in current PDB entries Yes, in about 100.0 % of entries
        # could also do the check manually (i.e. exists a chain with >= 2 nucleotides, but I wouldn't know for sure if they form a polymer)
        #       in the paper they allow only single nucleotides
        logger.warning(f'could not determine polymers of structure {s_id}, structure allowed, might contain RNA or DNA, though')

    # skip structure with too short chains
    if len(get_chains(s)) < 1:
        logger.info(f'skipping structure {s_id}: not enough chains')
        return False

    return True


def retrieve_structure_file_from_pdb(pdb_code: str) -> str:
    """ Download the structure file and return its path. If the file was already downloaded, use that.

    Downloads are in STRUCTURE_STORAGE_DIRECTORY. The caching is done by Bio.PDB.PDBList.

    :param pdb_code: pdb code
    :return: file name
    """

    return PDBList(
        pdb=Settings.STRUCTURE_STORAGE_DIRECTORY, verbose=logging.INFO >= logger.level
    ).retrieve_pdb_file(pdb_code, file_format='mmCif')
    # (mmCif is the default for file_format, but implicit = warning, this way no warning)
    # verbose: BioPython prints to stdout, determine verbose by root logger level (I could also have a settings variable
    # which also could be reset by a commandline argument)


def parse_structure(structure_file, structure_code=None):
    # bipythoní parser neni ideální, využívá legacy PDB fieldy (ATOM/HETATM) a auth_seq_id (s fallbackem na label_seq_id == problém při mapování z pdbe api), auth_asym_id. Pokud by např.
    # nebyl u HETATM auth_seq_id (nepovinný ale, in about 100.0 % of entries), spadlo by to
    # auth_seq_id u heteroatomů (v mmcifu mají všechny heteroatomy label_seq_id `.`) umožnuje identifikaci, do jaké molekuly atom patří (jinak by byl jen název sloučeniny)
    # ovšem ty auth_ položky nemusí být číslem, ale např. tento parser je převádí na int() -- může spadnout
    #    zde https://bioinformatics.stackexchange.com/a/11590 jsem se ale dočetl undocumented fakt,
    #    že ty písmenka nějak dali pryč a převedli je do insertion code sloupce
    # Anyway, I need label_seq_id of Residue objects not only for correct mapping of pdbe API results (see above), but also to align whole
    # polypeptide sequence (even unobserved residues) to form valid chain (apo-holo) pairs.

    # todo stejny jako input.get_structure, refactor

    mmcif_parser = CustomMMCIFParser()

    # remove warnings like: PDBConstructionWarning: WARNING: Chain C is discontinuous at line
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = mmcif_parser.get_structure(structure_code, structure_file)

    mmcif_dict_undocumented = mmcif_parser._mmcif_dict

    return structure, mmcif_parser.header, mmcif_dict_undocumented, BiopythonToMmcifResidueIds.create(mmcif_parser._mmcif_dict)


# def parse_and_filter_group_of_structures(files: Iterable[Any]) -> Iterable[Structure]:
#     """
#
#     :param files: Iterable of filenames or file-like objects of mmCIF files
#     :return: structures passing the filter
#     """
#     structure__header__mmcif_dict = map(parse_structure, files)
#
#     structures = (s for s, s_header, mmcif_dict in filter(structure_meets_our_criteria, structure__header__mmcif_dict))
#
#     return structures


def get_structure_filenames_in_dir(root_directory):
    for root, dirs, files in os.walk(root_directory):
        for name in files:
            yield os.path.join(root, name)


def print_random_groups(n):
    """ samples uniprot groups

    Provides test input, pick suitable groups and paste them below into dict `groups`.
    Run in interactive console, so that all_uniprot_groups is loaded only once (takes long), while sampling can be run multiple times """

    from apo_holo_structure_stats.input.uniprot_groups import get_basic_uniprot_groups
    all_uniprot_groups = get_basic_uniprot_groups()

    import random
    groups_sample = dict(random.sample(list(all_uniprot_groups.items()), n))
    print(groups_sample)


def get_test_input(groups=None):
    """ returns comma-separated sequence of structure PDB codes

    groups variable is like output from get_basic_uniprot_groups, or print_random_groups """

    from collections import defaultdict

    if groups is None:
        groups = {
            'Q2YQQ9': defaultdict(list,
                         {'3mqd': ['A'],
                          '3lrf': ['A'],
                          '4jv3': ['A'],
                          '3u0e': ['A'],
                          '3u0f': ['A']}),
            'A0ZZH6': defaultdict(list,  # z těhlech vyhovujou jen 3 holo: 5mb2 5man 5m9x
                                  {'2gdu': ['A', 'B'],
                                   '6fme': ['A', 'B'],
                                   '1r7a': ['A', 'B'],
                                   '5mb2': ['B'],
                                   '5c8b': ['B'],
                                   '5man': ['B'],
                                   '2gdv': ['A', 'B'],
                                   '5m9x': ['B']}),
        }

    # remove chain ids, keep just pdb codes
    struct_code_groups = [list(uniprot_group.keys()) for uniprot_group in groups.values()]

    import itertools
    return ','.join(itertools.chain(*struct_code_groups))


def get_chains_metadata_for_structure(ordinal, filename: Path, input_chain_metadata: dict = None,
                                      chain_whitelist: Set[str] = None):
    if input_chain_metadata is not None and chain_whitelist is None:
        raise ValueError('To update existing metadata about chains, chains have to be whitelisted.')

    parsed, s_metadata = parse_mmcif(path=filename, with_extra=True)
    s = parsed.structure

    logger.info(f'processing {ordinal}-th structure {s.id} at {filename}')

    # get the first model (s[0]), (in x-ray structures there is probably a single model)
    s, residue_id_mappings = map(lambda s: s[0], (s, parsed.bio_to_mmcif_mappings))

    get_chains_analyzer = GetChains()
    is_holo_analyzer = IsHolo()

    chain_metadata = []
    if structure_meets_our_criteria(s, s_metadata.header, s_metadata.mmcif_dict, get_chains_analyzer):
        for chain in get_chains_analyzer(s):
            # skip chains not in `chain_whitelist`
            if chain_whitelist and chain.id not in chain_whitelist:
                continue

            mapping_for_chain = residue_id_mappings[chain.id]

            # skip sequences with micro-heterogeneity
            # I think there were some phosphorylated residues in these cases, with altloc for that residue.
            # 50% P and non-P. But it doesn't make sense to compare the structure then (with itself or across).
            # As only the single residue has different coords.
            if mapping_for_chain.entity_poly_id in parsed.poly_with_microhet:
                logger.info(f'skipping {s.get_parent().id} {chain.id}. Microheterogeneity in sequence')
                continue

            # get chain metadata
            sequence = list(parsed.poly_seqs[mapping_for_chain.entity_poly_id].values())
            is_holo = is_holo_analyzer(s, ChainResidues.from_chain(chain, mapping_for_chain))

            try:
                metadata = input_chain_metadata[chain.id] if input_chain_metadata else {}
            except KeyError:
                # could fail, if chain.id not in chain_metadata (user error)
                raise RuntimeError(f'Missing metadata for chain {s.get_parent().id},{chain.id} in chain_metadata')

            metadata.update({
                'pdb_code': s.get_parent().id,
                'path': str(filename),  # todo I don't use that anymore
                'chain_id': chain.id,
                'is_holo': is_holo,  # possibly add some details about the ligand or LBS (which?)
                'resolution': get_resolution(s_metadata.mmcif_dict),
                '_exptl.method': s_metadata.mmcif_dict['_exptl.method'][0],  #   # there are structures where there are
                # multiple _exptl.method, but save just the first..
                'sequence': sequence})
            chain_metadata.append(metadata)

    return chain_metadata


def find_structures(pdb_codes):
    return map(lambda c: find_or_download_structure(c, allow_download=False), pdb_codes)


def main():
    import argparse
    import sys

    parser = get_argument_parser()
    parser.add_argument('--workers', type=int, default=1, help='number of subprocesses')
    parser.add_argument('--download_threads', type=int, default=1, help='number of threads')
    parser.add_argument('--all_chains', default=False, action='store_true')
    parser.add_argument('--disallow_download', default=False, action='store_true')
    parser.add_argument('--disallow_download', default=False, action='store_true')

    # parser.add_argument('--pdb_dir', type=str, action='store_true',
    #                     help='now pdb_codes_or_directory is a path to a directory with mmcif files. Whole tree structure is inspected, all files are assumed to be mmcifs.')
    parser.add_argument('-i', '--input_type', default='json', choices=['json', 'pdb_dir', 'pdb_codes'], help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('input', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('output_file', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')

    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    assert args.input_type == 'json'  # todo temporary hack (so that contains uniprotkb_id metadata)
    # todo nemusí, nevadí, pak si to tam může dodat, pokud to bude make_pairs vyžadovat

    chain_whitelists = None

    # translate input into structure filenames
    if args.input_type == 'pdb_dir':
        directory = args.input

        if not os.path.isdir(directory):
            logger.error(f'Directory {directory} does not exist')
            sys.exit(1)

        structure_filenames = get_structure_filenames_in_dir(directory)
    elif args.input_type == 'pdb_codes':
        pdb_codes = args.input.strip().split(',')

        if not pdb_codes:
            logger.error('No pdb codes specified')
            sys.exit(1)

        if args.disallow_download:
            structure_filenames = find_structures(pdb_codes)
        else:
            structure_filenames = download_structures(pdb_codes, args.download_threads)
        # structure_filenames = (retrieve_structure_file_from_pdb(pdb_code) for pdb_code in pdb_codes)
    elif args.input_type == 'json':
        # todo to accept just list of pdb_codes, add column chain_id? And won't that break chain_whitelists (want empty set)
        chains = pd.read_json(args.input)
        # chains = chains.iloc[:100]  # todo test hack
        # chimeric - more UNP to one chain (or could be in-vivo chimeric?
        # e.g. 2bnq E or 2bnu B have multiplem UNPs mapped to themselves
        # skip them all, (using one unp does not make sense) Or put it them with both unps?
        # todo obviously doesn't work, if this file is run with batches...
        chains = chains.drop_duplicates(subset=['pdb_code', 'chain_id'], keep=False)
        chains_gb_pdb_code = chains.groupby('pdb_code')

        metadata_gb_structure = chains_gb_pdb_code.apply(lambda df: df.set_index('chain_id').to_dict(orient='index'))
        chain_whitelists = chains_gb_pdb_code['chain_id'].apply(lambda series: set(series.to_list()))
        pdb_codes = chains_gb_pdb_code.indices.keys()
        if args.disallow_download:
            structure_filenames = find_structures(pdb_codes)
        else:
            structure_filenames = download_structures(pdb_codes, args.download_threads)
    else:
        raise ValueError('Unknown input type argument')

    structure_filenames = list(structure_filenames)
    logger.info(f'total structures to process: {len(structure_filenames)}')

    # load and filter structures
    extra_args = [metadata_gb_structure]
    if not args.all_chains and chain_whitelists is not None:
        extra_args.append(chain_whitelists)

    chains_of_structures_that_passed = []

    # serial version
    if args.workers == 1:
        for prepared_args in zip(itertools.count(), structure_filenames, *extra_args):
            try:
                chain_metadata = get_chains_metadata_for_structure(*prepared_args)
                chains_of_structures_that_passed.extend(chain_metadata)
            except Exception:
                logger.exception(f'Exception when a task for a structure `{prepared_args[1]}` was executed.')

    else:
        # parallel version, but logging somehow doesn't work and there is no easy/out-of-the-box way to make it work
        # with multiprocessing
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            chain_metadata_futures = submit_tasks(
                executor, 40 * args.workers,
                get_chains_metadata_for_structure, itertools.count(), structure_filenames, *extra_args
            )

            # iterate over the futures and flatten the chain metadata
            # result of a single task is a list of chain metadata for each structure, flatten tasks results into a list of chains
            for struct_filename, chains_future in zip(structure_filenames, chain_metadata_futures):
                print('done')
                try:
                    chain_metadata = chains_future.result()
                except Exception:
                    logger.exception(f'Exception when a task for a structure `{struct_filename}` was executed.')
                    continue

                # flatten
                chains_of_structures_that_passed.extend(chain_metadata)

    # with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
    #
    #
    #     # map to list of futures, so I can handle exceptions (with exeutor.map whole iterator stops in that case)
    #     chain_metadata_futures = list(map(
    #         lambda *args: executor.submit(get_chains_metadata_for_structure, *args),
    #         itertools.count(), structure_filenames, *extra_args,
    #     ))

    # iterate over the futures and flatten the chain metadata
    # result of a single task is a list of chain metadata for each structure, flatten tasks results into a list of chains
    # chains_of_structures_that_passed = []
    # for struct_filename, chains_future in zip(structure_filenames, chain_metadata_futures):
    #     try:
    #         chain_metadata = chains_future.result()
    #     except Exception:
    #         logger.exception(f'Exception when a task for a structure `{struct_filename}` was executed.')
    #         continue
    #
    #     # flatten
    #     chains_of_structures_that_passed.extend(chain_metadata)

    with open(args.output_file, 'w') as f:
        json.dump(chains_of_structures_that_passed, f)


if __name__ == '__main__':
    main()
