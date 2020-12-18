#!/usr/bin/env python3
import itertools
import json
from typing import List

from apo_holo_structure_stats.core.analyses import RMSD
from apo_holo_structure_stats.core.base_analyses import Analyzer, SerializableCachedAnalyzer


# serializer - u resultu a argumentů? Nebo u Analyzeru
# u argumentů by se to hodilo -- nemusel bych se opakovat u Analyzeru, u resultu taky
# ale u argumentů musí být i názvy -- a to asi chci konfigurovatelný, ne hackovat pythoní názvy parametrů
    # analyzer::serializátor může klidně volat na argumentech serialize (buď funkce, nebo metoda -- to bych pak musel patchnout v PBD.StructureBuilder importy SMCRA tříd), mixovat to, no nevim
        # bude vracet dict


class AnalysisRunner:
    def run(self, analyzer: Analyzer, *args, **kwargs):
        # jenze nekdy chci Analyzer jako Serializable
        raise NotImplementedError



class JSONAnalysisSerializer(AnalysisRunner):
    def __init__(self, output_file):
        self.output_file = output_file

    def run(self, analyzer: SerializableCachedAnalyzer, *args, **kwargs):
        # tady se prave typ analyzer parametru lisi (v pythonu to asi bude fungovat (mro) - typ nezajimavy, v jinych jazycich by se to neprelozilo, protoze run by neimplementovalo AnalysisRunner abstraktni metodu, neshoda typu argumentu)
        result = analyzer(*args, **kwargs)
        o = analyzer.serialize(result, *args, **kwargs)
        json.dump(o, self.output_file)
        # Unlike pickle and marshal, JSON is not a framed protocol, so trying to serialize multiple objects with repeated calls to dump() using the same fp will result in an invalid JSON file.



class ApoHoloFeeder:
    analyzers: List[Analyzer]

    def __init__(self):
        self.analyzers = []

    def add_analyzers(self, *analyzers: [Analyzer, ...]) -> None:
        self.analyzers.extend(analyzers)

    def run(self, apo_codes, holo_codes, get_structure):
        for apo_code, holo_code in itertools.product(apo_codes, holo_codes):
            apo, holo = map(get_structure, (apo_code, holo_code))

            for analyzer in self.analyzers:
                # spust job

                result = analyzer(apo, holo)



def run_analyses_for_isoform_group(apo_codes: List[str], holo_codes: List[str], get_structure):
    # apo-holo analyses
    apo_holo_pairs = itertools.product(apo_codes, holo_codes)

    apo_holo_runner = ApoHoloRunner()

    rmsd_a = RMSD()

    apo_holo_runner.add_analyzers(rmsd_a, rmsd_a)

    analyzers = []

    for apo_code, holo_code in apo_holo_pairs:
        apo, holo = map(get_structure, (apo_code, holo_code))



        print(apo.id, holo.id)
        print(rmsd_if_seqs_same(apo, holo))

    # holo-holo analyses
    holo_holo_pairs = itertools.combinations(holo_codes, 2)

    for holo1_code, holo2_code in holo_holo_pairs:
        print(holo1.id, holo2.id)
        print(rmsd_if_seqs_same(holo1, holo2))


if __name__ == '__main__':
    # run for all isoforms default
    # optionally specify a single isoform

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--isoform', help='process only structures with main chain of that isoform')
    parser.add_argument('structures_json_file', help='list of structures {pdb_code: , path: , isoform_id: , is_holo: bool, ?main_chain_id: }')
    parser.add_argument('output_file', help='dumped results of analyses')
    args = parser.parse_args()

    with open(args.structures_json_file) as f:
        structures_info = json.load(f)

    if args.isoform is not None:
        structures_info = list(filter(lambda s: s['isoform_id'] == args.isoform, structures_info))


    with open(args.output_file, 'w') as f:
        # serializace analyz bude asi prubezna

        # groupby isoform
            # run_analyses_for_isoform_group(, get_structure = lambda: 1  )
        pass
