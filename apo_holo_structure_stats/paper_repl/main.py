import concurrent.futures
import gzip
import itertools
import json
import logging
import os
import time
import warnings
from datetime import datetime
from difflib import SequenceMatcher
from multiprocessing import Manager
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List, Dict, Iterable

import pandas as pd
from Bio.PDB import PDBList, MMCIFParser

from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.Structure import Structure

from apo_holo_structure_stats.core.analyses import DomainResidues
from apo_holo_structure_stats.core.dataclasses import ChainResidueData, DomainResidueData, ChainResidues, DomainResidueMapping
from apo_holo_structure_stats.core.analysesinstances import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.input.download import APIException, download_and_save_file
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds, ResidueId, BioResidueId
from apo_holo_structure_stats.pipeline.run_analyses import chain_to_polypeptide, aligner, AnalysisHandler, \
    JSONAnalysisSerializer, ConcurrentJSONAnalysisSerializer


def get_longest_common_polypeptide(
        apo_seq: Dict[int, str], holo_seq: Dict[int, str],  # in 3-letter codes
    ):
    """ Substring """

    # find longest common substring, works even with the 3-letter codes (sequence of strings), e.g. MET and FME will get treated differently
    apo_residue_codes = list(apo_seq.values())
    holo_residue_codes = list(holo_seq.values())
    i1, i2, length = SequenceMatcher(a=apo_residue_codes, b=holo_residue_codes, autojunk=False)\
        .find_longest_match(0, len(apo_seq), 0, len(holo_seq))

    logging.info(f'substring to original ratio: {length/min(len(apo_seq), len(holo_seq))}')

    if length < MIN_SUBSTRING_LENGTH_RATIO * min(len(apo_seq), len(holo_seq)):
        alignment = next(aligner.align(apo_seq, holo_seq))
        logging.info('Sequences differ, alignment:')
        logging.info(f'\n{alignment}')
        logging.warning(f'does not meet the threshold for substring length')

    # crop polypeptides to longest common substring
    apo_common_seq = dict(itertools.islice(apo_seq.items(), i1, i1+length))
    holo_common_seq = dict(itertools.islice(holo_seq.items(), i2, i2+length))

    # return the sequences (values will be same, but not the keys, label_seq_ids, they might be offset, or depending on its definition
    # (which I find ambiguous), in a more complicated relation)
    return apo_common_seq, holo_common_seq


def get_observed_residues(
        chain1: Chain, c1_label_seq_ids: Iterable[int], c1_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
        chain2: Chain, c2_label_seq_ids: Iterable[int], c2_residue_mapping: BiopythonToMmcifResidueIds.Mapping):

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

        c1_residues.append(chain1[r1_bio_id])
        c2_residues.append(chain2[r2_bio_id])

        c1_residue_ids.append(r1_seq_id)
        c2_residue_ids.append(r2_seq_id)

    return c1_residues, c1_residue_ids, c2_residues, c2_residue_ids


MIN_SUBSTRING_LENGTH_RATIO = 0.90

def compare_chains(chain1: Chain, chain2: Chain,
                   c1_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c2_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c1_seq: Dict[int, str], c2_seq: Dict[int, str],  # in 3-letter codes
                   comparators__residues_param: List[Analyzer],
                   comparators__residue_ids_param: List[Analyzer],
                   comparators__domains__residues_param: List[Analyzer],
                   comparators__domains__residue_ids_param: List[Analyzer],
                   comparators__2domains__residues_param: List[Analyzer],
                       serializer_or_analysis_handler: AnalysisHandler,
                   domains_info: list,
                   ) -> None:
    """ Runs comparisons between two chains. E.g. one ligand-free (apo) and another ligand-bound (holo).
    :param chain1: A Bio.PDB Chain, obtained as a part of BioPython Structure object as usual
    :param chain2: A corresponding chain (same sequence), typically from a different PDB structure. See chain1.

    :param c1_residue_mapping:
    :param apo_poly_seqs:
    """
    s1_pdb_code = chain1.get_parent().get_parent().id
    s2_pdb_code = chain2.get_parent().get_parent().id

    logging.info(f'running analyses for ({s1_pdb_code}, {s2_pdb_code}) pair...')
    #
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     pp1 = chain_to_polypeptide(chain1)
    #     pp2 = chain_to_polypeptide(chain2)

    # c1_seq, c2_seq todo, is the order in atom_site loop guaranteed? If not, I should sort the dict by label_seq_id
    # also todo, is label_seq_id sequential, that is one-by-one always +1?
    # todo assert entity_poly_seq have no gaps (always +1), they say they're sequential, I think they mean exactly this

    # crop polypeptides to longest common substring
    c1_common_seq, c2_common_seq = get_longest_common_polypeptide(c1_seq, c2_seq)
    c1_label_seq_ids = list(c1_common_seq.keys())
    c2_label_seq_ids = list(c2_common_seq.keys())

    label_seq_id_offset = c2_label_seq_ids[0] - c1_label_seq_ids[0]

    # up to this point, we have residue ids of the protein sequence in the experiment. This also includes unobserved
    # residues, but those we will exclude from our analysis as their positions weren't determined
    c1_residues, c1_label_seq_ids, c2_residues, c2_label_seq_ids = get_observed_residues(
        chain1,
        c1_label_seq_ids,
        c1_residue_mapping,
        chain2,
        c2_label_seq_ids,
        c2_residue_mapping,
    )

    c1_residues = ChainResidues(c1_residues, s1_pdb_code, chain1.id)
    c2_residues = ChainResidues(c2_residues, s2_pdb_code, chain2.id)

    # todo trochu nesikovny
    c1_residue_ids = ChainResidueData[ResidueId]([ResidueId(label_seq_id, chain1.id) for label_seq_id in
                                                  c1_label_seq_ids], s1_pdb_code, chain1.id)
    c2_residue_ids = ChainResidueData[ResidueId]([ResidueId(label_seq_id, chain2.id) for label_seq_id in
                                                  c2_label_seq_ids], s2_pdb_code, chain2.id)

    # [done] tady nahradit pp pomocí apo_seq nějak
    # [done] v analyzerech (APIs) nahradit author_seq_id
    # todo tady matchovaní domén pomocí tohodle - zas mohu pouzit Sequence Matcher
    #   - ale spany, je to složitější -> zatím přeindexovat apo nebo holo do druhý...

    for a in comparators__residues_param:
        # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
        # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains

        serializer_or_analysis_handler.handle('chain2chain', a, a(c1_residues, c2_residues), c1_residues,
                                              c2_residues)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?
        # because what I would like is to output the analysis with objects identifiers, and then output the objects, what they contain (e.g. domain size?)


    for c in comparators__residue_ids_param:
        serializer_or_analysis_handler.handle('chain2chain', c, c(c1_residue_ids, c2_residue_ids), c1_residue_ids,
                                              c2_residue_ids)

    # domain-level analyses

    # get domains (set of auth_seq_id), sort them by domain id and hope they will correspond to each other
    # or could map corresponding domains by choosing the ones that have the most overlap?
    try:
        c1_domains = sorted(filter(lambda d: d.chain_id == chain1.id, get_domains(s1_pdb_code)), key=lambda d: d.domain_id)
        c2_domains = sorted(filter(lambda d: d.chain_id == chain2.id, get_domains(s2_pdb_code)), key=lambda d: d.domain_id)
        # todo zaznamenat total počet domén (pro obě struktury), zapsat do jinýho jsonu třeba
        domains_info.append(
            {'type': 'total_domains_found', 'result': len(c1_domains), 'pdb_code': s1_pdb_code, 'chain_id': chain1.id})
        domains_info.append(
            {'type': 'total_domains_found', 'result': len(c2_domains), 'pdb_code': s2_pdb_code, 'chain_id': chain2.id})
        # todo  spany domén, hlavně

    except APIException as e:
        if e.__cause__ and '404' in str(e.__cause__):
            logging.warning(f'{s1_pdb_code} {s2_pdb_code} no domains found, skip the domain-level analysis')
            return  # no domains found, skip the domain-level analysis
        raise


    # assert len(c1_domains) == len(c2_domains) # not always true, as expected, but now OK

    # SequenceMatcher on domain resiudes
    c1_domains__residues = []
    c2_domains__residues = []

    for c1_d in c1_domains:  # or c2_domains:
        # first remap first domain to second (or in future use longest common substrings, but not trivial since domains can be composed of multiple segments)
        # offset nemusí být všude stejný
        c1_domain_mapped_to_c2 = DomainResidueMapping.from_domain_on_another_chain(c1_d, chain2.id, label_seq_id_offset)

        # todo proc chain.get_parent??
        c1_d_residues = DomainResidues.from_domain(c1_d, chain1.get_parent(), c1_residue_mapping,
                                                   lambda id: id not in c1_label_seq_ids)
        c2_d_residues = DomainResidues.from_domain(c1_domain_mapped_to_c2, chain2.get_parent(), c2_residue_mapping,
                                                   lambda id: id not in c2_label_seq_ids)

        if not c1_d_residues or not c2_d_residues:
            # the domain is not within the processed LCS of both chains (empty intersection with chain residues)
            logging.warning(f'domain {c1_d.domain_id} is not within the processed LCS of both chains (empty '
                            f'intersection with '
                            f'chain residues)')
            continue

        # todo reflect domain cropping in the object id (domain id) somehow?
        c1_domains__residues.append(DomainResidues(c1_d_residues.data, c1_d_residues.structure_id, c1_d_residues.chain_id, c1_d_residues.domain_id))
        c2_domains__residues.append(DomainResidues(c2_d_residues.data, c2_d_residues.structure_id, c2_d_residues.chain_id, c2_d_residues.domain_id))

    # todo zaznamenat počet domén jdoucích do analýz
    domains_info.append({'type': 'analyzed_domain_count', 'result': len(c1_domains__residues), 'pdb_code': s1_pdb_code, 'chain_id': chain1.id})
    domains_info.append({'type': 'analyzed_domain_count', 'result': len(c2_domains__residues), 'pdb_code': s2_pdb_code, 'chain_id': chain2.id})

    # todo to tam taky neni v argumentech, ale harcoded.., to je ten muj fix...
    for chain_domains in (c1_domains__residues, c2_domains__residues):
        for d1, d2 in itertools.combinations(chain_domains, 2):
            serializer_or_analysis_handler.handle('2DA', get_interdomain_surface, get_interdomain_surface(d1, d2),
                                                  d1, d2)


    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        for a in comparators__domains__residues_param:
            serializer_or_analysis_handler.handle('domain2domain', a, a(d_chain1, d_chain2), d_chain1, d_chain2)

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

        for a in comparators__domains__residue_ids_param:
            serializer_or_analysis_handler.handle('domain2domain', a, a(d_chain1, d_chain2), d_chain1, d_chain2)

    # two-domain arrangements to two-domain arrangements
    for (d1_chain1, d1_chain2), (d2_chain1, d2_chain2) in itertools.combinations(zip(c1_domains__residues, c2_domains__residues), 2):
        # (in paper considered if of both apo and holo interdomain iface >= 200 A^2
        # if get_interdomain_surface(d1_chain1, d2_chain1) < 200 or get_interdomain_surface(d1_chain2, d2_chain2) < 200:
        #     continue

        for a in comparators__2domains__residues_param:
            serializer_or_analysis_handler.handle('chain2DA2chain2DA', a, a(d1_chain1, d2_chain1, d1_chain2,
                                                                            d2_chain2),
                                                  d1_chain1,
                                                  d2_chain1, d1_chain2, d2_chain2)

        d1d2_chain1 = d1_chain1 + d2_chain1
        d1d2_chain2 = d1_chain2 + d2_chain2
        serializer_or_analysis_handler.handle('chain2DA2chain2DA', get_rmsd, get_rmsd(d1d2_chain1, d1d2_chain2),
                                              d1d2_chain1,
                                              d1d2_chain2)  # todo hardcoded analysis

        # chain2DA2chain2DA nema stejny argumenty, asi v pohode, to je jenom pro level a moznost vybrat analyzu
        #   na danym levelu..

# def download_structure(pdb_code: str) -> str:
#     # logging.info(f'downloading structure {pdb_code}..')
#
#     filename = None
#
#     for i in range(1):
#         try:
#             filename = PDBList(pdb='pdb_structs', obsolete_pdb='pdb_structs', verbose=False).retrieve_pdb_file(pdb_code, file_format='mmCif')
#             break
#         except OSError as e:
#             if 'Too many connections' in str(e):
#                 print('too many connections')  # furt too many connections, i když stahuju přes 1 thread. Asi si to pamatuje dlouho..
#                 time.sleep(4)
#                 continue
#
#             # structure obsolete
#             try:
#                 logging.info(f'trying dl structure {pdb_code} as obsolete...')
#                 filename = PDBList(pdb='pdb_structs', obsolete_pdb='pdb_structs', verbose=False).retrieve_pdb_file(pdb_code, file_format='mmCif', obsolete=True)  # 1seo is obsolete, but we want to replicate the paper?
#                 logging.info('success')
#                 break
#             except OSError as e:
#                 logging.warning(f'{pdb_code} error: ' + str(e))
#
#     if filename is None:
#         logging.exception(f'downloading {pdb_code} failed: ')
#         raise
#     # else:
#     #     print(filename)
#     #     logging.info('no more too many connections, success')
#
#     # logging.info(f'finished downloading structure {pdb_code}')
#     return filename
#
#
#
# def get_structure(pdb_code: str) -> Structure:
#     filename = download_structure(pdb_code)
#     return MMCIFParser(QUIET=True).get_structure(pdb_code, filename)

def find_or_download_structure(pdb_code: str) -> str:
    """ Download structure file and return its path. Use existing file, if path already exists. """
    filename = f'{pdb_code}.cif.gz'
    url = f'https://files.rcsb.org/download/{filename}'

    local_path = Path('pdb_structs') / filename

    download_and_save_file(url, local_path)

    logging.info(f'finished downloading structure {pdb_code}, to {local_path}')

    return local_path


def get_structure(pdb_code: str):
    local_path = find_or_download_structure(pdb_code)

    with gzip.open(local_path, 'rt', newline='', encoding='utf-8') as text_file:
        mmcif_parser = MMCIFParser(QUIET=True)
        return mmcif_parser.get_structure(pdb_code, text_file), \
               BiopythonToMmcifResidueIds.create(mmcif_parser._mmcif_dict)  # reuse already parsed mmcifdict (albeit undocumented)


def get_chain_by_chain_code(model: Model, paper_chain_code: str) -> Chain:
    if paper_chain_code == '_':
        # seems that this is the code when there's only a single chain (probably named A?)
        chain = next(iter(model))

        if len(model) != 1 or chain.id != 'A':
            logging.warning(f'{model.get_parent().id} {paper_chain_code}, underscore is not what I thought. {chain.id}'
                            f', {len(model)}')
            long_chains = get_chains(model)
            if len(long_chains) > 1:
                logging.warning(f'model contains {len(long_chains)} chains with 50+ aa')
            return long_chains[0]  # sometime there is also chain B with ligands.

        return chain
    return model[paper_chain_code]


if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)

    def process_pair(s1_pdb_code: str, s2_pdb_code: str, s1_paper_chain_code: str, s2_paper_chain_code: str,
                     serializer, domains_info: list):

        # raise ValueError('hovna')
        logging.info(f'{s1_pdb_code}, {s2_pdb_code}')

        try:
            apo, (apo_residue_id_mappings, apo_poly_seqs) = get_structure(s1_pdb_code)
            holo, (holo_residue_id_mappings, holo_poly_seqs) = get_structure(s2_pdb_code)

            # get the first model (s[0]), (in x-ray structures there is probably a single model)
            apo, apo_residue_id_mappings = map(lambda s: s[0], (apo, apo_residue_id_mappings))
            holo, holo_residue_id_mappings = map(lambda s: s[0], (holo, holo_residue_id_mappings))

            apo_chain = get_chain_by_chain_code(apo, s1_paper_chain_code)
            holo_chain = get_chain_by_chain_code(holo, s2_paper_chain_code)

            apo_mapping = apo_residue_id_mappings[apo_chain.id]
            holo_mapping = holo_residue_id_mappings[holo_chain.id]

            compare_chains(apo_chain, holo_chain,
                           apo_mapping, holo_mapping,
                           apo_poly_seqs[apo_mapping.entity_poly_id], holo_poly_seqs[holo_mapping.entity_poly_id],
                           [get_rmsd],
                           [get_ss],
                           [get_rmsd],
                           [get_ss],
                           [get_hinge_angle],
                           serializer,
                           domains_info,
                           )
        except Exception as e:
            logging.exception('compare chains failed with: ')

    def run_apo_analyses():
        df = pd.read_csv('apo_holo.dat', delimiter=r'\s+', comment='#', header=None,
                         names=('apo', 'holo', 'domain_count', 'ligand_codes'), dtype={'domain_count': int})

        with Manager() as multiprocessing_manager:
            start_datetime = datetime.now()

            serializer = ConcurrentJSONAnalysisSerializer(f'output_apo_holo_{start_datetime.isoformat()}.json',
                                                          multiprocessing_manager)
            domains_info = multiprocessing_manager.list()

            found = False

            # with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
                for index, row in df.iterrows():
                    # if row.apo[:4] == '1cgj':
                    #     continue  # chymotrypsinogen x chymotrypsin + obojí má ligand... (to 'apo' má 53 aa inhibitor)

                    # if row.apo[:4] == '1ikp':
                    #     found = True
                    #
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
            with open(f'output_domains_info{start_datetime.isoformat()}.json', 'w') as f:
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

            apo = get_structure(row.apo[:4])
            holo = get_structure(row.holo[:4])

            apo_chain = get_chain_by_chain_code(apo, row.apo[4:])
            holo_chain = get_chain_by_chain_code(holo, row.holo[4:])

            try:
                compare_chains(apo_chain, holo_chain,
                               [get_rmsd],
                               [get_ss],
                               [get_rmsd],
                               [get_ss],
                               [get_hinge_angle],
                               serializer)
            except Exception as e:
                logging.exception('compare chains failed with: ')
                # raise

        pass


    run_apo_analyses()
    # run_holo_analyses()


