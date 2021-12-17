#!/usr/bin/env python3
import itertools
import json
import logging
import queue
from multiprocessing.managers import SyncManager
from typing import List, TypeVar, Generic, Iterable, Dict

from Bio.PDB import MMCIFParser, PPBuilder, is_aa
from Bio.PDB.Chain import Chain

from apo_holo_structure_stats import project_logger
from apo_holo_structure_stats.core.analyses import GetRMSD, GetMainChain, GetChains, CompareSecondaryStructure, \
    GetSecondaryStructureForStructure, GetDomainsForStructure, GetInterfaceBuriedArea, GetSASAForStructure, \
    GetCAlphaCoords, GetCentroid, GetCenteredCAlphaCoords, GetHingeAngle, GetRotationMatrix
from apo_holo_structure_stats.core.dataclasses import ChainResidueData, ChainResidues, DomainResidues, \
    DomainResidueMapping, DomainResidueData
from apo_holo_structure_stats.core.biopython_to_mmcif import ResidueId, BiopythonToMmcifResidueIds
from apo_holo_structure_stats.core.base_analyses import SerializableAnalyzer, Analyzer
from apo_holo_structure_stats.core.json_serialize import CustomJSONEncoder
from apo_holo_structure_stats.input.download import APIException, parse_mmcif
from apo_holo_structure_stats.pipeline.utils.log import add_loglevel_args

from apo_holo_structure_stats.core.analysesinstances import *


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


# debug aligner for preliminary sequence analysis (apo-holo/holo-holo), to see how they differ, if they differ
from Bio import Align
aligner = Align.PairwiseAligner(mode='global',
                                open_gap_score=-0.5,
                                extend_gap_score=-0.1,
                                end_gap_score=0,)  # match by default 1 and mismatch 0



# def chain_to_polypeptide(chain):
#     ppb = PPBuilder()
#     polypeptides = ppb.build_peptides(chain, aa_only=0)  # allow non-standard aas?
#
#     if len(polypeptides) != 1:
#         logging.info(f'parsed {len(polypeptides)} polypeptides from one chain, concatenating')
#
#         for pp in polypeptides[1:]:
#             polypeptides[0].extend(pp)
#
#     return polypeptides[0]

#
# def sequences_same(ch1, ch2):
#     pp1, pp2 = map(chain_to_polypeptide, (ch1, ch2))
#     seq1, seq2 = map(lambda pp: pp.get_sequence(), (pp1, pp2))
#
#     if seq1 != seq2:
#         # debug print to see how seqs differ
#
#         #     residues_mapping = # might be a fallback to the pdb-residue-to-uniprot-residue mapping api, depends, what we want
#         #           (also we have some mapping (segment to segment) when we do get_best_isoform, there can probably be mutations though)
#         #
#         #          if there are some extra residues on ends, we can truncate the sequences
#         #          maybe I wouldn't use the API, just do this alignment, if we wanted to compare structs with non-100%-identical sequences)
#
#         alignment = next(aligner.align(seq1, seq2))
#         logging.info('Sequences differ, alignment:')
#         logging.info(alignment)
#         return False
#
#     return True




def run_analyses_for_isoform_group(apo_codes: List[str], holo_codes: List[str], get_structure, serializer_or_analysis_handler: AnalysisHandler):
    # apo-holo analyses

    get_main_chain = GetMainChain((GetChains(),))

    get_c_alpha_coords = GetCAlphaCoords()
    get_centroid = GetCentroid((get_c_alpha_coords,))
    get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
    get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))

    get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
    get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

    ss_a = CompareSecondaryStructure((GetSecondaryStructureForStructure(),))
    interdomain_surface_a = GetInterfaceBuriedArea((GetSASAForStructure(),))

    comparators_of_apo_holo__residues_param = [get_rmsd, interdomain_surface_a]
    comparators_of_apo_holo__residue_ids_param = [ss_a]

    comparators_of_apo_holo_domains__residues_param = [get_rmsd, interdomain_surface_a]
    comparators_of_apo_holo_domains__residue_ids_param = [ss_a]

    comparators_of_apo_holo_2domains__residues_param = [get_hinge_angle]


    get_domains = GetDomainsForStructure()
    # nejrychlejsi /pro analyzu rmsd/ bude nacist CA z chainu/domen do np arraye. Pak porovnávat jen ty arraye, ale pak to nebude moc extensible pro jiny analyzy
    # napr interdomenovy povrch
    # tohle to zrychli víc nez nejaky caching centroidu (ten skoro vubec, rychly numpy, kdyz uz mam atomy v arrayi)
        # protoze se do pameti vejde vsechno a nemusi se pak znovu nacitat struktury (nebo nam bude stacit min pameti)
    # vubec bych si z tech analyz mozna mel vyzobnout na zacatku ty argumenty do mejch objektu

    # runner, co to umí spustit, jak chce, (plánuje multivláknově), podle těch kombinací dvojic třeba, jak má, je k ničemu, pokud ty argumenty budou plný, ,velký, objekty
    # protoze se zadela pamět jestě pred tim, nez se neco spusti. -> identifikatory chainů, domén (domény získám z apicka)

    # apo-holo analyses

    for apo_code, holo_code in itertools.product(apo_codes, holo_codes):
        logging.info(f'(not) running analyses for ({apo_code}, {holo_code}) apo-holo pair...')

        apo, holo = map(get_structure, (apo_code, holo_code))

        apo_main_chain, holo_main_chain = map(get_main_chain, (apo, holo))

        if not sequences_same(apo_main_chain, holo_main_chain):
            logging.info('skipping pair {apo_code}, {holo_code}. Main chain sequences differ.')
            continue

        apo_chain_residues, holo_chain_residues = map(
            lambda chain: ChainResidues([r for r in chain.get_residues() if is_aa(r)], chain.get_parent().get_parent().id, chain.id),
            (apo_main_chain, holo_main_chain)
        )

        for a in comparators_of_apo_holo__residues_param:
            # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
            # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains

            serializer_or_analysis_handler.handle('', a, a(apo_chain_residues, holo_chain_residues), apo_chain_residues,
                                                  holo_chain_residues)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?
            # because what I would like is to output the analysis with objects identifiers, and then output the objects, what they contain (e.g. domain size?)


        apo_chain_residue_ids, holo_chain_residue_ids = map(
            lambda chain_residues: ChainResidueData[ResidueId]([ResidueId.from_bio_residue(r) for r in chain_residues], chain_residues.structure_id, chain_residues.chain_id),
            (apo_chain_residues, holo_chain_residues)
        )

        for c in comparators_of_apo_holo__residue_ids_param:
            serializer_or_analysis_handler.handle('', c, c(apo_chain_residue_ids, holo_chain_residue_ids),
                                                  apo_chain_residue_ids, holo_chain_residue_ids)

        # domain analysis, asi by mely byt stejny (v sekvenci hopefully)
        #apo_domains = [] #apo.get_domains()  # opět může být cachované, tentokrát to bude malá response z apicka, obdobně SS
        # přesně! Tady je to nenačte, ty atomy. Pouze vrátí identifikátory a rozsah. Pokud bude někdo chtít, bude si moct to přeložit do BioPythoní entity, ten analyzer ale cachovaný nebude mezi fóry nebo vůbec
        #holo_domains = [] #holo.get_domains()
        # ještě se bude muset translatovat na array coordinates (to bude taky pomalý, ale nebude obrovský -- odhad
        # domena max 200, takze 200*3*8(double)= 4.8 kB= nic

        apo_domains = sorted(get_domains(apo_code), key=lambda d: d.domain_id,)
        holo_domains = sorted(get_domains(holo_code), key=lambda d: d.domain_id)


        # todo assert domains equal in sequence
        #   now equal in seqres? and correspond (only checks auth_seq_id/label_seq_id correspond, but that might be wrong - a different
        #   pdb structure might have a different numbering! TODO check sequence somehow (ideally without requiring to load the structure)
        for d_apo, d_holo in zip(apo_domains, holo_domains):
            assert len(d_apo) == len(d_holo) == sum(i == j for i,j in zip(d_apo, d_holo))

        apo_domains__residues = [DomainResidues.from_domain(d_apo, apo) for d_apo in apo_domains]
        holo_domains__residues = [DomainResidues.from_domain(d_holo, holo) for d_holo in holo_domains]

        for d_apo, d_holo in zip(apo_domains__residues, holo_domains__residues):

            for a in comparators_of_apo_holo_domains__residues_param:
                serializer_or_analysis_handler.handle('', a, a(d_apo, d_holo), d_apo, d_holo)

        for d_apo, d_holo in zip(apo_domains, holo_domains):
            d_apo = d_apo.to_set_of_residue_ids(apo_code)
            d_holo = d_holo.to_set_of_residue_ids(holo_code)

            for a in comparators_of_apo_holo_domains__residue_ids_param:
                serializer_or_analysis_handler.handle('', a, a(d_apo, d_holo), d_apo, d_holo)

        for (d1_apo, d1_holo), (d2_apo, d2_holo) in itertools.combinations(zip(apo_domains__residues,
                                                                               holo_domains__residues), 2):
            for a in comparators_of_apo_holo_2domains__residues_param:
                serializer_or_analysis_handler.handle('', a, a(d1_apo, d2_apo, d1_holo, d2_holo), d1_apo, d2_apo,
                                                      d1_holo, d2_holo)

            d1d2_apo = d1_apo + d2_apo
            d1d2_holo = d1_holo + d2_holo
            serializer_or_analysis_handler.handle('', get_rmsd, get_rmsd(d1d2_apo, d1d2_holo), d1d2_apo,
                                                  d1d2_apo)  # todo hardcoded analysis

    # holo-holo analyses

    h_h_struct_analyzers = [get_rmsd, ss_a]  # SS

    h_h_domain__analyzers = [get_rmsd]  # SS
    h_h_domain_pair_analyzers = [get_rmsd]  # rotation, screw axis, interdomain surface

    for holo1_code, holo2_code in itertools.combinations(holo_codes, 2):
        logging.info(f'running analyses for ({holo1_code}, {holo2_code}) holo-holo pair...')

        holo1, holo2 = map(get_structure, (holo1_code, holo2_code))

        # todo copy the preparation from apo-holo
        # for a in h_h_struct_analyzers:
        #     # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
        #     # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains
        #
        #     serializer_or_analysis_handler.handle(a, a(apo_chain_residues, holo_chain_residues), apo_chain_residues,
        #                                           holo_chain_residues)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?

        # domain analysis, asi by mely byt stejny (v sekvenci hopefully)
        holo1_domains = []#holo1.get_domains()  # vsechno může být již nacachované z apo-holo analýzy (otázka, todo jak jsou velké největší uniprot skupiny struktur, jestli se to vejde do paměti)
        holo2_domains = []#holo2.get_domains()

        corresponding_domains = list(zip(holo1_domains, holo2_domains))

        for d_holo1, d_holo2 in corresponding_domains:
            h_h_domain__analyzers

        for (d1_holo1, d1_holo2), (d2_holo1, d2_holo2) in itertools.combinations(corresponding_domains, 2):
            h_h_domain_pair_analyzers







# ADDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


def compare_chains(chain1: Chain, chain2: Chain,
                   c1_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c2_residue_mapping: BiopythonToMmcifResidueIds.Mapping,
                   c1_seq: Dict[int, str], c2_seq: Dict[int, str],  # in 3-letter codes
                   comparators__residues_param: List[Analyzer],  # todo tohle v pythonu 3.8 neni,
                   # ale spis bych mel mit duck-typing type annotations, protože to muze byt wrapply v cachich/a (de)serializerech..
                   # a nebude to inheritance ale composition (napr nebudu muset delat subclassy pro kazdy analyzer na to, abych tam pridal cache funkcionalitu...
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
    # also todo, is label_seq_id sequential, that is one-by-one always +1? - I use that contiguity assumption already...
    # - todo assert in code entity_poly_seq have no gaps (always +1), they say they're sequential, I think they mean exactly this
    #    - I could have done in filter structures, just to know for sure. If it were extensible already

    # crop polypeptides to longest common substring
    # todo - vzit z inputu (mam i1 a i2, staci ziskat offset a real label seq)
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

        for pdb_code, domains in ((s1_pdb_code, c1_domains), (s2_pdb_code, c2_domains)):
            for d in domains:
                domains_info.append(
                    {'type': 'full_domain',
                     'full_id': (pdb_code, d.chain_id, d.domain_id),
                     'pdb_code': pdb_code,
                     'chain_id': d.chain_id,
                     'domain_id': d.domain_id,
                     'spans': d.get_spans(),})


        # for d in c2_domains:
        #         domains_info.append(
        #     {'type': 'total_domains_found', 'result': len(c2_domains), 'pdb_code': s2_pdb_code, 'chain_id': chain2.id})
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
            logging.warning(f'domain {c1_d.domain_id} is not within the processed LCS of both chains (empty '
                            f'intersection with '
                            f'chain residues)')
            continue

        c1_domains__residues.append(DomainResidues(c1_d_residues.data, c1_d_residues.structure_id, c1_d_residues.chain_id, c1_d_residues.domain_id))
        c2_domains__residues.append(DomainResidues(c2_d_residues.data, c2_d_residues.structure_id, c2_d_residues.chain_id, c2_d_residues.domain_id))

    for residue_mapping, domains in ((c1_residue_mapping, c1_domains__residues),
                                     (c2_residue_mapping, c2_domains__residues)):
        for d in domains:
            domains_info.append(
                {'type': 'analyzed_domain',
                 'full_id': d.get_full_id(),
                 'pdb_code': d.structure_id,
                 'chain_id': d.chain_id,
                 'domain_id': d.domain_id,
                 'spans': d.get_spans(residue_mapping),
                 'spans_auth_seq_id': d.get_spans(residue_mapping, auth_seq_id=True),
                 })

    #
    # # todo zaznamenat počet domén jdoucích do analýz
    # domains_info.append({'type': 'analyzed_domain_count', 'result': len(c1_domains__residues), 'pdb_code': s1_pdb_code, 'chain_id': chain1.id})
    # domains_info.append({'type': 'analyzed_domain_count', 'result': len(c2_domains__residues), 'pdb_code': s2_pdb_code, 'chain_id': chain2.id})

    # todo to tam taky neni v argumentech, ale harcoded.., to je ten muj fix...
    # todo tohle totiž neni párový porovnání.., ale 'jednotkový'
    #  - stejně jako get domains, get_ss (nikoliv compare ss), vlastne i sequence atp
    #  - cachovat surface area teda nedava smysl, nacte se proste z predvypocitanyho, jako normalne
    #  - nebo, proste jenom tyhle structure-level veci ma smysl "cachovat" resp nepocitat tady, pro kazdej par, ale
    #  - nacitat z filu/unpicklovat - to asi ne, mít serialize/deserialize (stejne chci to mit jako citelny vystup). 4
    #  -  A pak to klidně všechno pro rychlost deserializovat do pameti...
    # no, tak to abych se těšil zas na json/pandas-merge hell.. Vsude merge.. Vsude dupe cols/delat index (ten pak ale nekdy zas potrebujes v cols...)
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
            serializer_or_analysis_handler.handle('chain2DA2chain2DA', a, a(d1_chain1, d2_chain1, d1_chain2, d2_chain2),
                                                  d1_chain1, d2_chain1, d1_chain2, d2_chain2)

        d1d2_chain1 = d1_chain1 + d2_chain1
        d1d2_chain2 = d1_chain2 + d2_chain2
        serializer_or_analysis_handler.handle('chain2DA2chain2DA', get_rmsd, get_rmsd(d1d2_chain1, d1d2_chain2),
                                              d1d2_chain1, d1d2_chain2)  # todo hardcoded analysis

        # chain2DA2chain2DA nema stejny argumenty, asi v pohode, to je jenom pro level a moznost vybrat analyzu
        #   na danym levelu..

#
# def get_chain_by_chain_code(model: Model, paper_chain_code: str) -> Chain:
#     if paper_chain_code == '_':
#         # seems that this is the code when there's only a single chain (probably named A?)
#         chain = next(iter(model))
#
#         if len(model) != 1 or chain.id != 'A':
#             logging.warning(f'{model.get_parent().id} {paper_chain_code}, underscore is not what I thought. {chain.id}'
#                             f', {len(model)}')
#             long_chains = get_chains(model)
#             if len(long_chains) > 1:
#                 logging.warning(f'model contains {len(long_chains)} chains with 50+ aa')
#             return long_chains[0]  # sometimes there is also chain B with ligands.
#
#         return chain
#     return model[paper_chain_code]


# todo je divny, ze se musi posilat ten Manager do kazdyho tasku a ne treba jen do toho ProcessPoolu
#    - tady asi ok, task trva celkem dlouho
#

# todo tohle uz normalne bezi (a pobezi) v jinym procesu, takze vsechny ty analyzery musim picklovat??
#  - nestaci je jenom definovat v globalu? - Nektery mozna, ale ty, co budou mít Manager (pdbe API analyzery)
#  - tak ty ne. Ty musim poslat, bohuzel, zejo
def process_pair(s1_pdb_code: str, s2_pdb_code: str, s1_paper_chain_code: str, s2_paper_chain_code: str,
                 serializer, domains_info: list):

    logging.info(f'{s1_pdb_code}, {s2_pdb_code}')

    try:
        apo, (apo_residue_id_mappings, apo_poly_seqs) = parse_mmcif(s1_pdb_code)
        holo, (holo_residue_id_mappings, holo_poly_seqs) = parse_mmcif(s2_pdb_code)
        # todo zkopirovat usage novyho API z filter structures..

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















def main():
    # runs for all isoforms by default
    # optionally specify a single isoform with --isoform

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--isoform', help='process only structures with main chain of that isoform')
    parser.add_argument('structures_json', help='list of structures {pdb_code: , path: , isoform_id: , is_holo: bool, ?main_chain_id: }')
    parser.add_argument('output_dir', help='dumped results of analyses')
    add_loglevel_args(parser)

    args = parser.parse_args()
    project_logger.setLevel(args.loglevel)
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    if args.isoform is not None:
        structures_info = list(filter(lambda s: s['isoform_id'] == args.isoform, structures_info))

    # with open(args.output_file, 'w') as f:
    #     kdyby serializace analyz byla prubezna (csv rows/triples stream)

    # groupby isoform
    key_isoform = lambda struct_info: struct_info['isoform_id']
    structures_info.sort(key=key_isoform)  # must be sorted for itertools.groupby

    # run analyses for each isoform group separately
    for isoform_id, isoform_group in itertools.groupby(structures_info, key=key_isoform):
        isoform_structure_info_dict = {s['pdb_code']: s for s in isoform_group}  # store info (path to files) for `structure_getter` lambda

        structure_getter = lambda pdb_code: MMCIFParser().get_structure(pdb_code, isoform_structure_info_dict[pdb_code]['path'])

        # divide pdb_codes into apo and holo lists
        apo_codes = []
        holo_codes = []

        for s in isoform_structure_info_dict.values():
            if s['is_holo']:
                holo_codes.append(s['pdb_code'])
            else:
                apo_codes.append(s['pdb_code'])

        # run analyses with a serializer/analysis handler
        serializer = JSONAnalysisSerializer(args.output_file)
        run_analyses_for_isoform_group(apo_codes, holo_codes, structure_getter, serializer)

        serializer.dump_data()


if __name__ == '__main__':
    main()


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
