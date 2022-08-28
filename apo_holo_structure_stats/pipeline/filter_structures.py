"""  Filters structures and fetches metadata using the parsed mmcif structures.

To modify the script functionality, you can inherit class StructureProcessor. """

import concurrent
import itertools
import json
import os
import warnings
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path
from typing import Set, List, Iterable

import pandas as pd
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import settings
from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analyses import GetChains, IsHolo
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds
from apo_holo_structure_stats.input.download import find_or_download_structure, parse_mmcif, MmcifParseResult
from apo_holo_structure_stats.pipeline.download_structures import download_structures
from apo_holo_structure_stats.pipeline.utils.log import get_argument_parser
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks, SemiBlockingQueueExecutor
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


def structure_meets_our_criteria(s, mmcif_dict, get_chains: GetChains):
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


def get_structure_filenames_in_dir(root_directory) -> Iterable[Path]:
    for root, dirs, files in os.walk(root_directory):
        for name in files:
            yield Path(root, name)


class StructureProcessor:
    """ Logic for structure filtering and fetching metadata from the parsed mmcif structures.

    To modify the script functionality, you can inherit this class and override some methods
    (retrieve_metadata_for_chain, keep_structure, ...) and set Settings.FILTER_STRUCTURES_CLASS to your descendant.
    """
    def __init__(self):
        """
        In case `args.workers` > 1, StructureProcessor.process_structure is run in multiple processes.
        You can use the __init__ to establish shared state with multiprocessing Manager if needed in your descendants
        of this class.
        """
        self.get_chains = GetChains()
        self.is_holo = IsHolo()  # this may not be needed in superclasses, but they still will probably call super __init__

    def keep_structure(self, parse_result: MmcifParseResult, mmcif_dict: MMCIF2Dict,
                       existing_chain_metadata: list = None, chain_whitelist: set = None) -> bool:
        return structure_meets_our_criteria(
            parse_result.structure[0],  # always choosing the first model
            mmcif_dict,
            self.get_chains,
        )

    def keep_chain(self, chain: Chain, id_mapping_chain: BiopythonToMmcifResidueIds.Mapping,
                   parse_result: MmcifParseResult, mmcif_dict: MMCIF2Dict) -> bool:
        # skip sequences with micro-heterogeneity
        # I think there were some phosphorylated residues in these cases, with altloc for that residue.
        # 50% P and non-P. But it doesn't make sense to compare the structure then (with itself or across).
        # As only the single residue has different coords.
        if id_mapping_chain.entity_poly_id in parse_result.poly_with_microhet:
            logger.info(f'skipping {chain.get_parent().get_parent().id} {chain.id}. Microheterogeneity in sequence')
            return False

        return True

    def retrieve_metadata_for_chain(self, filename: Path, parse_result: MmcifParseResult, mmcif_dict: MMCIF2Dict,
                                    chain: Chain, id_mapping_chain: BiopythonToMmcifResidueIds.Mapping,
                                    existing_metadata: dict) -> dict:
        sequence = list(parse_result.poly_seqs[id_mapping_chain.entity_poly_id].values())
        is_holo = self.is_holo(chain.get_parent(), ChainResidues.from_chain(chain, id_mapping_chain))

        metadata = existing_metadata.copy()
        metadata.update({
            'pdb_code': chain.get_parent().get_parent().id,
            'path': str(filename),  # todo I don't use that anymore, wouldn't work later in pipeline
            'chain_id': chain.id,
            'is_holo': is_holo,  # possibly add some details about the ligand or LBS (which?)
            'resolution': get_resolution(mmcif_dict),
            '_exptl.method': mmcif_dict['_exptl.method'][0],  # # there are structures where there are
            # multiple _exptl.method, but save just the first..
            'sequence': sequence})
        return metadata

    def retrieve_chain_metadata(self, filename: Path, parse_result: MmcifParseResult, mmcif_dict: MMCIF2Dict,
                                existing_chain_metadata: list = None, chain_whitelist: set = None) -> Iterable[dict]:
        s = parse_result.structure[0]   # always choosing the first model
        id_mappings_chain = parse_result.bio_to_mmcif_mappings[0]

        for chain in self.get_chains(s):
            # skip chains not in `chain_whitelist`
            if chain_whitelist and chain.id not in chain_whitelist:
                continue

            if not self.keep_chain(chain, id_mappings_chain[chain.id], parse_result, mmcif_dict):
                continue

            try:
                metadata = existing_chain_metadata[chain.id] if existing_chain_metadata else {}
            except KeyError:
                # could fail, if chain.id not in chain_metadata (user error)
                raise RuntimeError(f'Missing metadata for chain {s.get_parent().id},{chain.id} in chain_metadata')

            yield self.retrieve_metadata_for_chain(filename, parse_result, mmcif_dict, chain,
                                                   id_mappings_chain[chain.id], metadata)

    def process_structure(self, filename: Path, existing_chain_metadata: list = None,
                          chain_whitelist: set = None) -> List[dict]:
        if existing_chain_metadata is not None and chain_whitelist is None:
            raise ValueError('To update existing metadata about chains, chains have to be whitelisted.')

        # load structure
        parsed, s_metadata = parse_mmcif(path=filename, with_mmcif_dict=True)

        # filter out the structure or return retrieved metadata for each chain
        if self.keep_structure(parsed, s_metadata, existing_chain_metadata, chain_whitelist):
            return list(self.retrieve_chain_metadata(filename, parsed, s_metadata, existing_chain_metadata,
                                                     chain_whitelist))
        return []


def find_structures(pdb_codes):
    return map(lambda c: find_or_download_structure(c, allow_download=False), pdb_codes)


def main():
    import sys

    parser = get_argument_parser()
    parser.add_argument('--workers', type=int, default=1, help='number of subprocesses')
    parser.add_argument('--download_threads', type=int, default=1, help='number of threads')
    parser.add_argument('--disallow_download', default=False, action='store_true')

    parser.add_argument('-i', '--input_type', default='json', choices=['json', 'pdb_dir', 'pdb_codes'],
                        help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.'
                             'In that case, whole tree structure is inspected, all files are assumed to be mmcifs, or'
                             'gzipped mmcifs (filename ending with .gz)')
    parser.add_argument('input', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('output_file', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')

    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()

    assert args.input_type == 'json'  # todo temporary (so that contains uniprotkb_id metadata, and I don't use `path`
    # in metadata anywhere else. So -i pdb_dir wouldn't work - the structures couldn't be loaded in run_analyses then,
    # until I used `path` there.)

    extra_args = []

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
    elif args.input_type == 'json':
        chains = pd.read_json(args.input)
        # chimeric - more UNP to one chain (or could be in-vivo chimeric?
        # e.g. 2bnq E or 2bnu B have multiple UNPs mapped to themselves
        # skip them all, (using one unp does not make sense) Or put it them with both unps?
        # todo obviously doesn't work, if this file is run with batches...
        chains = chains.drop_duplicates(subset=['pdb_code', 'chain_id'], keep=False)
        chains_gb_pdb_code = chains.groupby('pdb_code')

        metadata_gb_structure = chains_gb_pdb_code.apply(lambda df: df.set_index('chain_id').to_dict(orient='index'))
        chain_whitelists = chains_gb_pdb_code['chain_id'].apply(lambda series: set(series.to_list()))

        extra_args = [metadata_gb_structure, chain_whitelists]

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
    chains_of_structures_that_passed = []

    # get the configurable class that filters the structures and returns metadata for chains
    sp: StructureProcessor = settings.load_class_from_str(Settings.FILTER_STRUCTURES_CLASS)()

    # serial version
    if args.workers == 1:
        for counter, filename, *extra in zip(itertools.count(), structure_filenames, *extra_args):
            logger.info(f'processing {counter + 1}-th structure at {filename}')
            try:
                chain_metadata = sp.process_structure(filename, *extra)
                chains_of_structures_that_passed.extend(chain_metadata)
            except Exception:
                logger.exception(f'Exception when a task for a structure `{filename}` was executed.')

    else:
        # parallel version, but logging somehow doesn't work and there is no easy/out-of-the-box way to make it work
        # with multiprocessing
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            chain_metadata_futures = submit_tasks(
                executor, 40 * args.workers,
                sp.process_structure, structure_filenames, *extra_args
            )

            # iterate over the futures and flatten the chain metadata
            # (result of a single task is a list of chain metadata for each structure,
            # flatten tasks results into a list of chains)
            for counter, filename, chains_future in zip(itertools.count(), structure_filenames, chain_metadata_futures):
                logger.info(f'processing {counter + 1}-th structure at {filename}')
                try:
                    chain_metadata = chains_future.result()
                except Exception:
                    logger.exception(f'Exception when a task for a structure `{filename}` was executed.')
                    continue

                # flatten
                chains_of_structures_that_passed.extend(chain_metadata)

    with open(args.output_file, 'w') as f:
        json.dump(chains_of_structures_that_passed, f)


if __name__ == '__main__':
    main()
