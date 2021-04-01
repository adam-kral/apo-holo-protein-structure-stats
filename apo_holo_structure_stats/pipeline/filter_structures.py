#!/usr/bin/env python3

import json
import os
import warnings
import logging
from typing import Iterable, List, Any

from Bio.PDB import PDBList, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Structure import Structure

from apo_holo_structure_stats.core.analyses import GetChains, GetMainChain
from apo_holo_structure_stats.pipeline.log import add_loglevel_args
from apo_holo_structure_stats.settings import STRUCTURE_DOWNLOAD_ROOT_DIRECTORY, MIN_STRUCTURE_RESOLUTION


def structure_meets_our_criteria(s, s_header, mmcif_dict, get_chains: GetChains):
    """ decides if structure meets criteria for resolution, single-chainedness, etc. """

    resolution = s_header['resolution']

    # skip low resolution
    if resolution and resolution > MIN_STRUCTURE_RESOLUTION:
        logging.info(f'skipping structure {s.id}: resolution ({resolution}) does not meet the limit of {MIN_STRUCTURE_RESOLUTION}')
        return False

    # skip non-xray. Done in the original paper. Try without it. Or then, specify in output dataset which experimental
    #   method (or one could integrate the info themselves)
    pass

    # skip DNA/RNA complexes
    try:
        if any(poly_type in mmcif_dict['_entity_poly.type'] for poly_type in ('polydeoxyribonucleotide', 'polyribonucleotide')):
            logging.info(f'skipping structure {s.id}: no interest in structures with RNA or DNA chains')
            return False
    except KeyError:
        # _entity_poly.type:: Required in PDB entries no; Used in current PDB entries Yes, in about 100.0 % of entries
        # could also do the check manually (i.e. exists a chain with >= 2 nucleotides, but I wouldn't know for sure if they form a polymer)
        #       in the paper they allow only single nucleotides
        logging.warning(f'could not determine polymers of structure {s.id}, structure allowed, might contain RNA or DNA, though')

    # skip non-single-chain structure (too short chain or too many chains)
    if len(get_chains(s)) < 1:
        logging.info(f'skipping structure {s.id}: not a enough chains')
        return False

    return True


class CustomMMCIFParser(MMCIFParser):
    """ Adapted BioPython code, just to get the pdb code (_entry.id) from the mmcif file """
    def get_structure(self, file, structure_id=None):
        """ Parses file contents and returns Structure object.

        Note that parameter order is different to the BioPython's implementation (reversed, as structure_id is optional).

        :param file: a file-like object or a file name
        :param structure_id: if not specified, taken from mmcif (`_entry.id`)
        :return: Bio.PDB.Structure
        """

        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            self._mmcif_dict = MMCIF2Dict(file)

            # begin change
            if structure_id is None:
                structure_id = self._mmcif_dict['_entry.id'][0].lower()
            # end change

            self._build_structure(structure_id)
            self._structure_builder.set_header(self._get_header())

        return self._structure_builder.get_structure()


def retrieve_structure_file_from_pdb(pdb_code: str) -> str:
    """ Download the structure file and return its path. If the file was already downloaded, use that.

    Downloads are in STRUCTURE_DOWNLOAD_ROOT_DIRECTORY. The caching is done by Bio.PDB.PDBList.

    :param pdb_code: pdb code
    :return: file name
    """

    return PDBList(pdb=STRUCTURE_DOWNLOAD_ROOT_DIRECTORY, verbose=logging.INFO >= logging.root.level).retrieve_pdb_file(pdb_code, file_format='mmCif')
    # (mmCif is the default for file_format, but implicit = warning, this way no warning)
    # verbose: BioPython prints to stdout, determine verbose by root logger level (I could also have a settings variable
    # which also could be reset by a commandline argument)


def parse_structure(structure_file, structure_code=None):
    # bipythoní parser neni ideální, využívá legacy PDB fieldy (ATOM/HETATM) a auth_seq_id (s fallbackem na label_seq_id == problém při mapování z pdbe api), auth_asym_id. Pokud by např.
    # nebyl u HETATM auth_seq_id (nepovinný ale, in about 100.0 % of entries), spadlo by to
    #  auth_asym_id, in about 93.8 % of entries -- zde to může spadnout
    # auth_seq_id u heteroatomů (v mmcifu mají všechny heteroatomy label_seq_id `.`) umožnuje identifikaci, do jaké molekuly atom patří (jinak by byl jen název sloučeniny)
    # ovšem ty auth_ položky nemusí být číslem, ale např. tento parser je převádí na int() -- může spadnout
    mmcif_parser = CustomMMCIFParser()

    # remove warnings like: PDBConstructionWarning: WARNING: Chain C is discontinuous at line
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = mmcif_parser.get_structure(structure_file, structure_code)

    mmcif_dict_undocumented = mmcif_parser._mmcif_dict

    return structure, mmcif_parser.header, mmcif_dict_undocumented


def parse_and_filter_group_of_structures(files: Iterable[Any]) -> Iterable[Structure]:
    """

    :param files: Iterable of filenames or file-like objects of mmCIF files
    :return: structures passing the filter
    """
    structure__header__mmcif_dict = map(parse_structure, files)

    structures = (s for s, s_header, mmcif_dict in filter(structure_meets_our_criteria, structure__header__mmcif_dict))

    return structures


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

# todo logging level ArgumentParser default subclass?

if __name__ == '__main__':
    # chce, aby se tomu mohly dodat fily přes directory
    # nebo se tomu dá comma-delimited list of pdb_codes

    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--is-directory', action='store_true',
                        help='now pdb_codes_or_directory is a path to a directory with mmcif files. Whole tree structure is inspected, all files are assumed to be mmcifs.')

    parser.add_argument('pdb_codes_or_directory', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('output_file', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')
    add_loglevel_args(parser)

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    # translate input into structure filenames
    if args.is_directory:
        directory = args.pdb_codes_or_directory

        if not os.path.isdir(directory):
            logging.error(f'Directory {directory} does not exist')
            sys.exit(1)

        structure_filenames = get_structure_filenames_in_dir(directory)
    else:
        pdb_codes = args.pdb_codes_or_directory.strip().split(',')

        if not pdb_codes:
            logging.error('No pdb codes specified')
            sys.exit(1)

        structure_filenames = (retrieve_structure_file_from_pdb(pdb_code) for pdb_code in pdb_codes)

    structure_filenames = list(structure_filenames)
    logging.info(f'total structures to process: {len(structure_filenames)}')

    # load and filter structures
    structures_passed = []

    for i, filename in enumerate(structure_filenames):
        s, s_header, mmcif_dict = parse_structure(filename)
        logging.info(f'processing {i}-th structure {s.id} at {filename}')

        get_chains_analyzer = GetChains()

        if structure_meets_our_criteria(s, s_header, mmcif_dict, get_chains_analyzer):
            structures_passed.append({'pdb_code': s.id, 'path': filename, 'main_chain_id': GetMainChain((get_chains_analyzer,))(s).id})  # maybe full path name?

            # path so we know where the structure file is in the next pipeline steps,
            # main_chain_id for isoform step (so we don't need to load the structure again just to know the main chain id)

    with open(args.output_file, 'w') as f:
        json.dump(structures_passed, f)
