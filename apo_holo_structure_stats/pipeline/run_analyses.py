#!/usr/bin/env python3
import itertools
import json
from typing import List, TypeVar, Generic

from Bio.PDB import MMCIFParser

from apo_holo_structure_stats.core.analyses import RMSD, GetMainChain, GetChains
from apo_holo_structure_stats.core.base_analyses import Analyzer, SerializableCachedAnalyzer, SerializableAnalyzer


TAnalyzer = TypeVar('TAnalyzer')


class AnalysisHandler(Generic[TAnalyzer]):
    def handle(self, analyzer: TAnalyzer, result, *args, **kwargs):
        raise NotImplementedError


class JSONAnalysisSerializer(AnalysisHandler[SerializableAnalyzer]):
    # json na to neni uplne vhodny, neda se dumpovat nebo nacitat inkrementalně, jako napr. csvcko -- je treba mit vzdy v pameti celou reprezentaci (vsechny vysledky, nez se to dumpne)
    def __init__(self, output_file_name):
        self.output_file_name = output_file_name
        self.data = []

    def handle(self, analyzer: SerializableAnalyzer, result, *args, **kwargs):
        """ serializes analysis args and results into a json"""

        o = analyzer.serialize(result, *args, **kwargs)
        self.data.append(o)
          # nejde dumpovat to asi nejde takhle inkrementálně. Musí se pak dumpnout asi nastřádaný všechno
        # Unlike pickle and marshal, JSON is not a framed protocol, so trying to serialize multiple objects with repeated calls to dump() using the same fp will result in an invalid JSON file.

    def dump_data(self):
        with open(self.output_file_name, 'w') as f:
            json.dump(self.data, f)


# class Feeder:
#     analyzers: List[Analyzer]
#     analysis_handler: AnalysisHandler  # muze se zmenit na List, kdyby bylo treba (output ve vice formatech). Pak by bylo rozumny posilat ale result a args a spoustet to zde
#
#     def __init__(self, analysis_handler: AnalysisHandler):
#         self.analyzers = []
#         self.analysis_handler = analysis_handler
#
#     def add_analyzers(self, *analyzers: [Analyzer, ...]) -> None:
#         self.analyzers.extend(analyzers)
#
#     def run_batch(self, *args, **kwargs):
#         raise NotImplementedError
#
#     def run_analyzers_with_args(self, *args, **kwargs):
#         for analyzer in self.analyzers:
#             result = analyzer(*args, **kwargs)
#
#             self.analysis_handler.handle(analyzer, result, *args, **kwargs)
#
#
# class ApoHoloPairFeeder(Feeder):
#     """ Feeds pairs into analyzers/  indirectly via analysishandlers """
#
#     def run_batch(self, apo_codes, holo_codes, get_structure):
#         # neni finalni, protoze tohle nemuze plug and play vyuzit parallelfeeder. Protože to getuje structure
#         for apo_code, holo_code in itertools.product(apo_codes, holo_codes):
#             apo, holo = map(get_structure, (apo_code, holo_code))
#
#             apo_domains = apo.get_domains()
#             holo_domains = holo.get_domains()
#             # asi by mely byt stejny (v sekvenci hopefully)
#
#             corresponding_domains = list(zip(apo_domains, holo_domains))
#
#             pairs =
#
#             # da se rict, ze zde jsou pouze twin struct analyzers (asi bych to mohl dat jako TwinStructAnalyzer?)
#             self.run_analyzers_with_args(apo, holo)


# vyhodit feeders, jsou k ničemu
# stejne bude jen jedna pipeline
# pokud bude multithread, bude se pouzivat ta


def run_analyses_for_isoform_group(apo_codes: List[str], holo_codes: List[str], get_structure, serializer_or_analysis_handler: AnalysisHandler):
    # apo-holo analyses

    rmsd_a = RMSD((GetMainChain((GetChains(),)),))

    a_h_struct_analyzers = [rmsd_a]  # SS
    a_h_domain_analyzers = [rmsd_a]  # SS
    a_h_domain_pair_analyzers = [rmsd_a]  # taky rmsd (jejich graf rmsd vs bending), rotation (vyuzije rmsd, tady to zrovna asi smysl dává), screw axis (posunutí), interdomain surface
        # tady se vyuzije caching? rotační matrix na zarovnání první domény

    # nejrychlejsi /pro analyzu rmsd/ bude nacist CA z chainu/domen do np arraye. Pak porovnávat jen ty arraye, ale pak to nebude moc extensible pro jiny analyzy
    # napr interdomenovy povrch
    # tohle to zrychli víc nez nejaky caching centroidu (ten skoro vubec, rychly numpy, kdyz uz mam atomy v arrayi)
        # protoze se do pameti vejde vsechno a nemusi se pak znovu nacitat struktury (nebo nam bude stacit min pameti)
    # vubec bych si z tech analyz mozna mel vyzobnout na zacatku ty argumenty do mejch objektu

    # runner, co to umí spustit, jak chce, (plánuje multivláknově), podle těch kombinací dvojic třeba, jak má, je k ničemu, pokud ty argumenty budou plný, ,velký, objekty
    # protoze se zadela pamět jestě pred tim, nez se neco spusti. -> identifikatory chainů, domén (domény získám z apicka)


    # apo-holo analyses

    for apo_code, holo_code in itertools.product(apo_codes, holo_codes):
        apo, holo = map(get_structure, (apo_code, holo_code))

        for a in a_h_struct_analyzers:
            # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
            # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains

            serializer_or_analysis_handler.handle(a, a(apo, holo), apo_code, holo_code)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?
            # because what I would like is to output the analysis with objects identifiers, and then output the objects, what they contain (e.g. domain size?)


    # ---------- odtud níže zatim nefunkční stub až do konce funkce

        # domain analysis, asi by mely byt stejny (v sekvenci hopefully)
        apo_domains = [] #apo.get_domains()  # opět může být cachované, tentokrát to bude malá response z apicka, obdobně SS
        holo_domains = [] #holo.get_domains()
        # ještě se bude muset translatovat na array coordinates (to bude taky pomalý, ale nebude obrovský -- odhad
        # domena max 200, takze 200*3*8(double)= 4.8 kB= nic

        corresponding_domains = list(zip(apo_domains, holo_domains))

        for d_apo, d_holo in corresponding_domains:
            a_h_domain_analyzers

        for (d1_apo, d1_holo), (d2_apo, d2_holo) in itertools.combinations(corresponding_domains, 2):
            a_h_domain_pair_analyzers

        # # da se rict, ze zde jsou pouze twin struct analyzers (asi bych to mohl dat jako TwinStructAnalyzer?)
        # self.run_analyzers_with_args(apo, holo)

    # holo-holo analyses

    h_h_struct_analyzers = [rmsd_a]  # SS

    h_h_domain__analyzers = [rmsd_a]  # SS
    h_h_domain_pair_analyzers = [rmsd_a]  # rotation, screw axis, interdomain surface

    for holo1_code, holo2_code in itertools.combinations(holo_codes, 2):
        holo1, holo2 = map(get_structure, (holo1_code, holo2_code))

        h_h_struct_analyzers

        # domain analysis, asi by mely byt stejny (v sekvenci hopefully)
        holo1_domains = []#holo1.get_domains()  # vsechno může být již nacachované z apo-holo analýzy (otázka, todo jak jsou velké největší uniprot skupiny struktur, jestli se to vejde do paměti)
        holo2_domains = []#holo2.get_domains()

        corresponding_domains = list(zip(holo1_domains, holo2_domains))

        for d_holo1, d_holo2 in corresponding_domains:
            h_h_domain__analyzers

        for (d1_holo1, d1_holo2), (d2_holo1, d2_holo2) in itertools.combinations(corresponding_domains, 2):
            h_h_domain_pair_analyzers


if __name__ == '__main__':
    # runs for all isoforms by default
    # optionally specify a single isoform with --isoform

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--isoform', help='process only structures with main chain of that isoform')
    parser.add_argument('structures_json', help='list of structures {pdb_code: , path: , isoform_id: , is_holo: bool, ?main_chain_id: }')
    parser.add_argument('output_file', help='dumped results of analyses')
    args = parser.parse_args()

    with open(args.structures_json) as f:
        structures_info = json.load(f)

    if args.isoform is not None:
        structures_info = list(filter(lambda s: s['isoform_id'] == args.isoform, structures_info))

    # with open(args.output_file, 'w') as f:
    #     kdyby serializace analyz byla prubezna (csv rows/triples stream)

    # groupby isoform
    key_isoform = lambda struct_info: struct_info['isoform_id']
    structures_info.sort(key=key_isoform)  # must be sorted for itertools.groupby

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
