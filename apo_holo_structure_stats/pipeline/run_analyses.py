import concurrent
import itertools
import json
import logging
import queue
import sys
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import Manager
from multiprocessing.managers import SyncManager
from pathlib import Path
from typing import List, TypeVar, Generic, Iterable, Dict, Any

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PPBuilder, is_aa
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analyses import GetRMSD, GetMainChain, GetChains, CompareSecondaryStructure, \
    GetSecondaryStructureForStructure, GetDomainsForStructure, GetInterfaceBuriedArea, GetSASAForStructure, \
    GetCAlphaCoords, GetCentroid, GetCenteredCAlphaCoords, GetHingeAngle, GetRotationMatrix, AnalysisException, \
    MissingDataException
from apo_holo_structure_stats.core.dataclasses import ChainResidueData, ChainResidues, DomainResidues, \
    DomainResidueMapping, DomainResidueData
from apo_holo_structure_stats.core.biopython_to_mmcif import ResidueId, BiopythonToMmcifResidueIds
from apo_holo_structure_stats.core.base_analyses import SerializableAnalyzer, Analyzer
from apo_holo_structure_stats.core.json_serialize import CustomJSONEncoder
from apo_holo_structure_stats.input.download import APIException, parse_mmcif
from apo_holo_structure_stats.pipeline.make_pairs_lcs import fn_wrapper_unpack_args, LCSResult, \
    pairs_without_mismatches, load_pairs_json
from apo_holo_structure_stats.pipeline.run_analyses_settings import configure_pipeline
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args, get_argument_parser

# from apo_holo_structure_stats.core.analysesinstances import *
from apo_holo_structure_stats.pipeline.utils.task_queue import submit_tasks

logger = logging.getLogger(__name__)

TAnalyzer = TypeVar('TAnalyzer')


class AnalysisHandler(Generic[TAnalyzer]):
    def handle(self, level_tag, analyzer: TAnalyzer, result, *args, **kwargs):
        raise NotImplementedError


class JSONAnalysisSerializer(AnalysisHandler[SerializableAnalyzer]):
    # json na to neni uplne vhodny, neda se dumpovat nebo nacitat inkrementalně, jako napr. csvcko -- je treba mit vzdy v pameti celou reprezentaci (vsechny vysledky, nez se to dumpne)
    def __init__(self, output_file_name):
        self.output_file_name = output_file_name
        self.data = []

    def handle(self, level_tag, analyzer: SerializableAnalyzer, result, *args, **kwargs):
        """ serializes analysis args and results into a json
        :param level_tag:
        """

        o = analyzer.serialize(result, *args, **kwargs)
        o['level_tag'] = level_tag
        self.data.append(o)
          # nejde dumpovat to asi nejde takhle inkrementálně. Musí se pak dumpnout asi nastřádaný všechno
        # Unlike pickle and marshal, JSON is not a framed protocol, so trying to serialize multiple objects with repeated calls to dump() using the same fp will result in an invalid JSON file.

    def dump_data(self):
        with open(self.output_file_name, 'w') as f:
            json.dump(self.data, f, cls=CustomJSONEncoder)


class ConcurrentJSONAnalysisSerializer(AnalysisHandler[SerializableAnalyzer]):
    """ Works with multiprocessing """
    # json na to neni uplne vhodny, neda se dumpovat nebo nacitat inkrementalně, jako napr. csvcko -- je treba mit vzdy v pameti celou reprezentaci (vsechny vysledky, nez se to dumpne)
    def __init__(self, output_file_name, multiprocessing_manager: SyncManager):
        self.output_file_name = output_file_name
        # self.queue = Queue()  # or I could use Manager.list(), but python docs say that queue should be preferred to
        self.queue = multiprocessing_manager.Queue()
        # sharing state (managers with lists),
        # I don't need to consume the queue concurrently, only once, at the end, so I use only the putting functionality

        # STACKOVERFLOW:
        # multiprocessing.Pool already has a shared result-queue, there is no need to additionally involve a
        # Manager.Queue. Manager.Queue is a queue.Queue (multithreading-queue) under the hood, located on a separate
        # server-process and exposed via proxies. This adds additional overhead compared to Pool's internal queue.

    def handle(self, level_tag, analyzer: SerializableAnalyzer, result, *args, **kwargs):
        """ serializes analysis args and results into a json
        :param level_tag:
        """

        o = analyzer.serialize(result, *args, **kwargs)
        o['level_tag'] = level_tag
        self.queue.put(o)
          # nejde dumpovat to asi nejde takhle inkrementálně. Musí se pak dumpnout asi nastřádaný všechno
        # Unlike pickle and marshal, JSON is not a framed protocol, so trying to serialize multiple objects with repeated calls to dump() using the same fp will result in an invalid JSON file.

    def dump_data(self):
        data = []
        try:
            while 1:
                data.append(self.queue.get_nowait())
        except queue.Empty:
            pass

        with open(str(self.output_file_name), 'w') as f:
            json.dump(data, f, cls=CustomJSONEncoder)


# process-local variable, to be set in each worker process with ProcessPoolExecutor
# main process sets it too
class PLocal:
    pass
# could be a dataclass
plocal = PLocal()


def assert_label_seq_id_contiguous(seq: List[int]):
    """Sanity check for the format of mmcif (spec not totally clear in this, therefore sanity check) """
    for i in range(len(seq) - 2):
        assert seq[i] + 1 == seq[i+1]


def get_seqs_range_and_offset(seq1: Dict[int, str], seq2: Dict[int, str], lcs_result: LCSResult):
    for seq in (seq1, seq2):
        try:
            assert_label_seq_id_contiguous(list(seq.keys()))
        except AssertionError:
            logger.warning('ERROR, label_seq list is not contiguous as expected. What now?')
            try:
                assert_label_seq_id_contiguous(sorted(seq.keys()))
            except AssertionError:
                logger.warning('ERROR, even if sorted, label_seq list is not contiguous as expected. What now?')

    def get_seq_range(seq, lcs_start):
        seq_ids = list(seq.keys())
        start = seq_ids[lcs_start]
        return range(start, start + lcs_result.length + 1)

    seq1_range = get_seq_range(seq1, lcs_result.i1)
    seq2_range = get_seq_range(seq2, lcs_result.i2)
    seq2_label_seq_id_offset = seq2_range[0] - seq1_range[0]
    return seq1_range, seq2_range, seq2_label_seq_id_offset


def compare_chains(chain1: Chain, chain2: Chain,
                   c1_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c2_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c1_seq: Dict[int, str], c2_seq: Dict[int, str],  # in 3-letter codes
                   lcs_result: LCSResult) -> None:
    """ Runs comparisons between two chains. E.g. one ligand-free (apo) and another ligand-bound (holo).
    :param chain1: A Bio.PDB Chain, obtained as a part of BioPython Structure object as usual
    :param chain2: A corresponding chain (same sequence), typically from a different PDB structure. See chain1.

    :param c1_residue_mapping:
    """
    s1_pdb_code = chain1.get_parent().get_parent().id
    s2_pdb_code = chain2.get_parent().get_parent().id

    logger.info(f'running analyses for ({s1_pdb_code}, {s2_pdb_code}) pair...')

    # c1_seq, c2_seq todo, is the order in atom_site loop guaranteed? If not, I should sort the dict by label_seq_id
    # also is label_seq_id sequential, that is one-by-one always +1? - I use that contiguity assumption already...
    # - the following asserts entity_poly_seq have no gaps (always +1), they say they're sequential, I think they mean exactly this
    #    - I could have done in filter structures, just to know for sure
    c1_seq_range, c2_seq_range, label_seq_id_offset = get_seqs_range_and_offset(c1_seq, c2_seq, lcs_result)

    # up to this point, we have residue ids of the protein sequence in the experiment. This also includes unobserved
    # residues, but those we will exclude from our analysis as their positions weren't determined
    c1_residues, c1_label_seq_ids, c2_residues, c2_label_seq_ids = get_observed_residues(
        chain1,
        list(c1_seq_range),
        c1_residue_mapping,
        chain2,
        list(c2_seq_range),
        c2_residue_mapping,
    )

    c1_residues = ChainResidues(c1_residues, s1_pdb_code, chain1.id)
    c2_residues = ChainResidues(c2_residues, s2_pdb_code, chain2.id)

    # todo trochu nesikovny
    c1_residue_ids = ChainResidueData[ResidueId]([ResidueId(label_seq_id, chain1.id) for label_seq_id in
                                                  c1_label_seq_ids], s1_pdb_code, chain1.id)
    c2_residue_ids = ChainResidueData[ResidueId]([ResidueId(label_seq_id, chain2.id) for label_seq_id in
                                                  c2_label_seq_ids], s2_pdb_code, chain2.id)

    def run_analysis(level_tag, analysis, *args):
        try:
            plocal.serializer_or_analysis_handler.handle(level_tag, analysis, analysis(*args), *args)
        except AnalysisException:
            logger.exception('Caught exception while computing analysis, all others will be run normally')

    # chain-level analyses
    for a in plocal.comparators__residues_param:
        # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
        # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains

        run_analysis('chain2chain', a, c1_residues, c2_residues)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?
        # because what I would like is to output the analysis with objects identifiers, and then output the objects, what they contain (e.g. domain size?)


    for a in plocal.comparators__residue_ids_param:
        run_analysis('chain2chain', a, c1_residue_ids, c2_residue_ids)

    # todo save the residue ids which were compared (already have this info in the analyzed_domain object below) but
    #   when the chain does not have domains available from API (low percentage though)

    # domain-level analyses

    # get domains (set of auth_seq_id), sort them by domain id and hope they will correspond to each other
    # or could map corresponding domains by choosing the ones that have the most overlap?
    try:
        c1_domains = sorted(filter(lambda d: d.chain_id == chain1.id, plocal.get_domains(s1_pdb_code)), key=lambda d: d.domain_id)
        c2_domains = sorted(filter(lambda d: d.chain_id == chain2.id, plocal.get_domains(s2_pdb_code)), key=lambda d: d.domain_id)

        # # todo tohle muzu dat uplne jinam.. treba do 1struct remote analyses, jestli to dá do ramky celej json vubec?
        # no ale ja musim nejak vyhodnotit ty analyzy v tom jupyter notebooku, jak to vůbec nactu do ramky úplně všechno vlastně??
        # asi nenačtu... Takže to musim zrušit z jupytera, nebo si udělám interactive job a pustim jupytera na metacentru s ramkou 48 gb treba a pripojim se na nej
        # jde tam vubec delat server?

        for pdb_code, domains in ((s1_pdb_code, c1_domains), (s2_pdb_code, c2_domains)):
            key = (plocal.get_domains.get_name(), pdb_code)

            if key not in plocal.one_struct_analyses_done_set:  # todo, tady jsem to zapomnel nastavit, takze to vlastne nepouzivam. Jsou tam dupes, ale je to jedno v podstate..
                for d in domains:
                    plocal.domains_info.append(
                        {'type': 'full_domain',
                         'full_id': (pdb_code, d.chain_id, d.domain_id),
                         'pdb_code': pdb_code,
                         'chain_id': d.chain_id,
                         'domain_id': d.domain_id,
                         'spans': d.get_spans(),})

    except APIException as e:
        if e.__cause__ and '404' in str(e.__cause__):
            logger.warning(f'{s1_pdb_code} {s2_pdb_code} no domains found, skip the domain-level analysis')
            return  # no domains found, skip the domain-level analysis
        raise
    except MissingDataException:
        logger.warning(f'{s1_pdb_code} {s2_pdb_code} no domains found, skip the domain-level analysis')
        return  # no domains found, skip the domain-level analysis


    # assert len(c1_domains) == len(c2_domains) # not always true, as expected, but now OK

    c1_domains__residues = []
    c2_domains__residues = []

    # todo zatim to necham, ale uz mam i to defakto domain lcs
    for c1_d in c1_domains:  # or c2_domains:
        # first remap first domain to second (or in future use longest common substrings, but not trivial since domains can be composed of multiple segments)
        # offset nemusí být všude stejný
        c1_domain_mapped_to_c2 = DomainResidueMapping.from_domain_on_another_chain(c1_d, chain2.id, label_seq_id_offset)

        # todo proc chain.get_parent?? Asi abych chain nemusel specifikovat (ale ted pracuju jenom s nima..)
        c1_d_residues = DomainResidues.from_domain(c1_d, chain1.get_parent(), c1_residue_mapping,
                                                   lambda id: id not in c1_label_seq_ids)
        c2_d_residues = DomainResidues.from_domain(c1_domain_mapped_to_c2, chain2.get_parent(), c2_residue_mapping,
                                                   lambda id: id not in c2_label_seq_ids)

        if not c1_d_residues or not c2_d_residues:
            # the domain is not within the processed LCS of both chains (empty intersection with chain residues)
            logger.warning(f'domain {c1_d.domain_id} is not within the processed LCS of both chains (empty '
                            f'intersection with '
                            f'chain residues)')
            continue

        c1_domains__residues.append(DomainResidues(c1_d_residues.data, c1_d_residues.structure_id, c1_d_residues.chain_id, c1_d_residues.domain_id))
        c2_domains__residues.append(DomainResidues(c2_d_residues.data, c2_d_residues.structure_id, c2_d_residues.chain_id, c2_d_residues.domain_id))

    for residue_mapping, domains in ((c1_residue_mapping, c1_domains__residues),
                                     (c2_residue_mapping, c2_domains__residues)):
        for d in domains:
            # samozrejme blbost to ukladat pokazdy, kdyz je to paru s necim klidne nekolikrat...
            # EDIT Ne uplne, je to croply na observed residues z obou párů (a lcs)
            # ale musel bych upravit idcko a dat tam nejak ten pár...
            # dupes pres joby muzu mazat pres ['pair_id', 'full_id']
            plocal.domains_info.append(
                {'type': 'analyzed_domain',
                 'pair_id': (c1_residues.get_full_id(), c2_residues.get_full_id()),  # domain cropping depends on the paired chains
                 'full_id': d.get_full_id(),
                 'pdb_code': d.structure_id,
                 'chain_id': d.chain_id,
                 'domain_id': d.domain_id,
                 'spans': d.get_spans(residue_mapping),
                 'spans_auth_seq_id': d.get_spans(residue_mapping, auth_seq_id=True),
                 })

    #  - cachovat surface area teda nedava smysl, nacte se proste z predvypocitanyho, jako normalne
    #  - nebo, proste jenom tyhle structure-level veci ma smysl "cachovat" resp nepocitat tady, pro kazdej par, ale
    #  -  A pak to klidně všechno pro rychlost deserializovat do pameti...
    # no, tak to abych se těšil zas na json/pandas-merge hell.. Vsude merge.. Vsude dupe cols/delat index (ten pak ale nekdy zas potrebujes v cols...)

    for chain_domains in (c1_domains__residues, c2_domains__residues):
        for d1, d2 in itertools.combinations(chain_domains, 2):
            # key = (plocal.get_interdomain_surface.get_name(),) + (d1.get_full_id(), d2.get_full_id())

            pair_id = (c1_residues.get_full_id(), c2_residues.get_full_id())
            # hack nepouzivam local fn run_analysis, chci tam i pair_id
            plocal.serializer_or_analysis_handler.handle('2DA', plocal.get_interdomain_surface, plocal.get_interdomain_surface(d1, d2),
                                                  d1, d2, pair_id)


    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        for a in plocal.comparators__domains__residues_param:
            run_analysis('domain2domain', a, d_chain1, d_chain2)

    # todo vyres ty divny idcka
    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        # Convert DomainResidues to DomainResidueData[ResidueId]
        # asi zas přes mapping... lepší by to bylo, kdyby byl implicitně schovaný třeba na to biopython residue (
        # jinak by to nešlo moc ani, leda mit CustomResidue s fieldama bioresidue a label_seq_id, to je ale celkem
        # naprd, nebo ne? Nefungovalo by to s chainem, ale to stejně nikde nepoužívám...
        d_chain1 = DomainResidueData[ResidueId]([ResidueId.from_bio_residue(r, c1_residue_mapping) for r in d_chain1],
                                                d_chain1.structure_id, d_chain1.chain_id, d_chain1.domain_id)
        d_chain2 = DomainResidueData[ResidueId]([ResidueId.from_bio_residue(r, c2_residue_mapping) for r in d_chain2],
                                                d_chain2.structure_id, d_chain2.chain_id, d_chain2.domain_id)

        for a in plocal.comparators__domains__residue_ids_param:
            run_analysis('domain2domain', a, d_chain1, d_chain2)

    # two-domain arrangements to two-domain arrangements
    for (d1_chain1, d1_chain2), (d2_chain1, d2_chain2) in itertools.combinations(zip(c1_domains__residues, c2_domains__residues), 2):
        # (in paper considered if of both apo and holo interdomain iface >= 200 A^2
        # if get_interdomain_surface(d1_chain1, d2_chain1) < 200 or get_interdomain_surface(d1_chain2, d2_chain2) < 200:
        #     continue

        for a in plocal.comparators__2DA__residues_param:
            run_analysis('chain2DA2chain2DA', a, d1_chain1, d2_chain1, d1_chain2, d2_chain2)


        d1d2_chain1 = d1_chain1 + d2_chain1
        d1d2_chain2 = d1_chain2 + d2_chain2
        run_analysis('chain2DA2chain2DA', plocal.get_rmsd, d1d2_chain1, d1d2_chain2)  # todo hardcoded analysis

        # chain2DA2chain2DA nema stejny argumenty, asi v pohode, to je jenom pro level a moznost vybrat analyzu
        #   na danym levelu..


def process_pair(s1_pdb_code: str, s2_pdb_code: str, s1_chain_code: str, s2_chain_code: str,
                 lcs_result: LCSResult, compare_chains_fn=compare_chains):
    logger.info(f'process_pair {s1_pdb_code}, {s2_pdb_code}')

    try:
        apo_parsed = plocal.parse_mmcif(s1_pdb_code, allow_download=True)
        holo_parsed = plocal.parse_mmcif(s2_pdb_code, allow_download=True)

        apo = apo_parsed.structure
        apo_residue_id_mappings = apo_parsed.bio_to_mmcif_mappings
        apo_poly_seqs = apo_parsed.poly_seqs

        holo = holo_parsed.structure
        holo_residue_id_mappings = holo_parsed.bio_to_mmcif_mappings
        holo_poly_seqs = holo_parsed.poly_seqs

        # get the first model (s[0]), (in x-ray structures there is probably a single model)
        apo, apo_residue_id_mappings = map(lambda s: s[0], (apo, apo_residue_id_mappings))
        holo, holo_residue_id_mappings = map(lambda s: s[0], (holo, holo_residue_id_mappings))

        apo_chain = apo[s1_chain_code]
        holo_chain = holo[s2_chain_code]

        apo_mapping = apo_residue_id_mappings[apo_chain.id]
        holo_mapping = holo_residue_id_mappings[holo_chain.id]

        compare_chains_fn(apo_chain, holo_chain,
                       apo_mapping, holo_mapping,
                       apo_poly_seqs[apo_mapping.entity_poly_id], holo_poly_seqs[holo_mapping.entity_poly_id],
                       lcs_result,
                       # *analyzers,
                       # get_domains,
                       # get_rmsd,
                       # get_interdomain_surface,
                       # serializer,
                       # domains_info,
                       # one_struct_analyses_done_set,
                       )
    except Exception as e:
        logger.exception('compare chains failed with: ')


# todo instance analýz, taky poslat cely dostany listy - ajaj, mám vubec ty instance v jinejch procesech?
#  Ne, musi musim je tam proste poslat z master procesu...
#   - later reply: stejné instance protože jsem chtěl mít to cachovaní, aby to fungovalo s tim manager.List atp.

# MP initializer
# https://stackoverflow.com/questions/10117073/how-to-use-initializer-to-set-up-my-multiprocess-pool

def run_analyses_serial(pairs, process_pair_fn, worker_initializer=None, worker_init_args=()):
    # pairs = pairs.iloc[:20]  # todo testing hack
    # print(pairs)
    if worker_initializer:
        worker_initializer(*worker_init_args)

    start_datetime = datetime.now()

    for i, row in enumerate(pairs.itertuples()):
        process_pair_fn(row.pdb_code_apo, row.pdb_code_holo,
               row.chain_id_apo, row.chain_id_holo,
               row.lcs_result)
        logger.info(f'done {i + 1} / {len(pairs)}')

    logger.info(f'run_analyses duration: {timedelta(seconds=int((datetime.now() - start_datetime).seconds))}\n'
                f'started at: {start_datetime.isoformat()}, ended at: {datetime.now().isoformat()}')


def run_analyses_multiprocess(pairs, process_pair_fn, num_workers, worker_initializer=None, worker_init_args=None):
    """ Allows to run `process_pair_fn` on pairs (in multiprocess environment).

    Reads the json in the expected format (output of `make_pairs_lcs.py`). Hides some of the multiprocessing complexity.

    If `process_pair_fn` should access and modify some shared state between processes, you might want to use params
    `worker_initializer` and `worker_init_args` to set up the shared state. (Possibly with the SyncManager derived
     objects -- shared lists, dicts, etc.) Then, this function would be called inside a `with multiprocessing.Manager()`
     block.

    :param pairs:
    :param process_pair_fn: picklable callable, taking 5 parameters: pdb_code_apo, _holo, chain_id_apo, _holo, lcs_result
    :param num_workers:
    :param worker_initializer: callable (probably) needs to be picklable, for more details look up this parameter
        in ProcessPoolExecutor constructor
    :param worker_init_args: arguments passed to worker_initializer in each worker, again picklable,  for more details
        look up `initargs` in ProcessPoolExecutor constructor
    :return:
    """

    start_datetime = datetime.now()

    # with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:    # todo testing
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=worker_initializer,
            initargs=worker_init_args) as executor:
        def get_args():
            for row in pairs.itertuples():
                yield (row.pdb_code_apo, row.pdb_code_holo,
                       row.chain_id_apo, row.chain_id_holo,
                       row.lcs_result,)

        fn = partial(fn_wrapper_unpack_args, process_pair_fn)
        futures = submit_tasks(executor, 40 * num_workers, fn, get_args())
        # wait for all futures to complete
        for i, f in enumerate(futures):
            f.result()
            # log progress
            logger.info(f'done {i + 1} / {len(pairs)}')

        logger.info(f'run_analyses duration: {timedelta(seconds=int((datetime.now() - start_datetime).seconds))}\n'
                    f'started at: {start_datetime.isoformat()}, ended at: {datetime.now().isoformat()}')


def worker_initializer(analyzer_namespace, serializer_or_analysis_handler, domains_info, one_struct_analyses_done_set):
    attrs = locals()
    del attrs['analyzer_namespace']
    attrs.update(analyzer_namespace)

    for attr_name, value in attrs.items():
        setattr(plocal, attr_name, value)


parser = get_argument_parser(description='Compares chains given the pairs.',
                             epilog='Currently not extensible. Users can write their own script, '
                                    'similar to this one, or redefine `configure_pipeline` function for smaller changes.')
parser.add_argument('--workers', default=1, type=int, help='Workers to use. Now only 1 is supported. '
                                                           'Currently, for "multiprocessing" just run the script '
                                                           'multiple times with chunks of data. That is actually _the_ '
                                                           'solution when run on a cluster with a job scheduler.')
parser.add_argument('--opt_input_dir', type=Path, default=Path())
# parser.add_argument('chains_json', help='list of structures {pdb_code: , path: , isoform_id: , is_holo: bool, ?main_chain_id: }')
parser.add_argument('pairs_json', help='Output of ah-make-pairs, list of JSON objects with fields: '
                                       'pdb_code_apo, pdb_code_holo, chain_id_apo, chain_id_holo, lcs_result, '
                                       'where lcs_result is serialized `make_pairs_lcs.LCSResult`')


def main():
    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)
    logging.basicConfig()

    assert args.workers == 1, 'todo: do multiprocessing differently (chunks) or use just one worker – ' \
                              'thats even better when run on a cluster'

    potential_pairs = load_pairs_json(args.pairs_json)
    if potential_pairs.empty:
        logger.warning('Input json contains no records.')
        sys.exit(0)
    pairs = pairs_without_mismatches(potential_pairs)

    # don't run analyses for each uniprot group separately, as creating a process pool carries an overhead
    # median pairs per group is 6
    # but could reset the caches? No need, all are LRU..
    start_datetime = datetime.now()
    analyses_output_fpath = Path(f'output_apo_holo_{start_datetime.isoformat()}.json')
    domains_info_fpath = Path(f'output_domains_info_{start_datetime.isoformat()}.json')

    with Manager() as multiprocessing_manager:
        # get analyzers as configured
        # p = configure_pipeline(multiprocessing_manager)
        analyses_namespace = configure_pipeline(args.opt_input_dir)

        serializer = ConcurrentJSONAnalysisSerializer(analyses_output_fpath, multiprocessing_manager)
        # domains_info = multiprocessing_manager.list()
        # one_struct_analyses_done_set = multiprocessing_manager.dict()
        domains_info = []
        one_struct_analyses_done_set = {}
        worker_initializer_args = (analyses_namespace, serializer, domains_info, one_struct_analyses_done_set)

        # run_analyses_multiprocess(pairs, process_pair, args.workers, worker_initializer,
        #                           worker_initializer_args)
        run_analyses_serial(pairs, process_pair, worker_initializer,
                                  worker_initializer_args)

        serializer.dump_data()
        with domains_info_fpath.open('w') as f:
            json.dump(list(domains_info), f)

# todo proc to bezi tak dlouho
#   - velky struktury viru treba 6v1z
#   - maj hodne paru (totiz hodne stejnejch chainu..)
#   - hlavne trvaj ty struktury desne dlouho nacist..


if __name__ == '__main__':
    main()


def c_alpha_observed(residue: Residue):
    return 'CA' in residue


def get_observed_residues(
        chain1: Chain, c1_label_seq_ids: Iterable[int], c1_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
        chain2: Chain, c2_label_seq_ids: Iterable[int], c2_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
        require_c_alphas_observed=True):
    """ Returns residues observed in both chains (and their label_seq_ids).

    Observed = there are coordinates for residue's atom (any) in the mmcif ATOM_SITE category.
    If `require_c_alphas_observed` is True, the residues need to contain coordinates for their c_alpha atom.
    """
    c1_residues = []
    c2_residues = []
    c1_residue_ids = []
    c2_residue_ids = []

    for r1_seq_id, r2_seq_id in zip(c1_label_seq_ids, c2_label_seq_ids):
        try:
            r1_bio_id = c1_residue_mapping.to_bio(r1_seq_id)
            r2_bio_id = c2_residue_mapping.to_bio(r2_seq_id)
        except KeyError:
            # a residue unobserved (wasn't in atom list) -> skip the whole pair
            continue

        r1 = chain1[r1_bio_id]
        r2 = chain2[r2_bio_id]

        if require_c_alphas_observed and not (c_alpha_observed(r1) and c_alpha_observed(r2)):
            continue

        c1_residues.append(r1)
        c2_residues.append(r2)

        c1_residue_ids.append(r1_seq_id)
        c2_residue_ids.append(r2_seq_id)

    return c1_residues, c1_residue_ids, c2_residues, c2_residue_ids
