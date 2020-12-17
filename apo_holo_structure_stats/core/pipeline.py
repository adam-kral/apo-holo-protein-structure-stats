import itertools
from collections import defaultdict, namedtuple
import logging
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser, is_aa, NeighborSearch, PDBList


from apo_holo_structure_stats.input.download import get_structure_stringio, get_best_isoform_for_chains
from apo_holo_structure_stats.input.uniprot_groups import get_basic_uniprot_groups




def get_chains(struct):
    return filter(lambda chain: sum(is_aa(residue) for residue in chain) >= 50, struct.get_chains())



def get_main_chain(struct):
    return next(get_chains(struct))  # todo cache get_main_chain


def get_isoform(struct):
    isoform_mapping = get_best_isoform_for_chains(struct.id)  # relying in pdb id in parsed structure object, not-so-good design?

    return isoform_mapping[get_main_chain(struct).id][0].uniprot_id  # getting the first segment's uniprot_id (could be a chimera, but probably won't)


# definitions for analyses
from Bio.PDB.Polypeptide import PPBuilder
import rmsd

def chain_to_polypeptide(chain):
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(chain, aa_only=0)  # allow non-standard aas?

    if len(polypeptides) != 1:
        print('warning ', len(polypeptides), ' polypeptides from one chain, extending first pp')

        for pp in polypeptides[1:]:
            polypeptides[0].extend(pp)

    return polypeptides[0]


# debug aligner for preliminary sequence analysis (apo-holo/holo-holo), to see how they differ, if they differ
from Bio import Align
aligner = Align.PairwiseAligner(mode='global',
                                open_gap_score=-0.5,
                                extend_gap_score=-0.1,
                                end_gap_score=0,)  # match by default 1 and mismatch 0


def get_c_alpha_coords(polypeptide):
    return np.array([res['CA'].get_coord() for res in polypeptide])


def rmsd_if_seqs_same(struct1, struct2):
    pp1, pp2 = map(lambda struct: chain_to_polypeptide(get_main_chain(struct)), (struct1, struct2))
    seq1, seq2 = map(lambda pp: pp.get_sequence(), (pp1, pp2))

    if seq1 != seq2:
        # debug print to see how seqs differ

        #     residues_mapping = # might be a fallback to the pdb-residue-to-uniprot-residue mapping api, depends, what we want
        #           (also we have some mapping (segment to segment) when we do get_best_isoform, there can probably be mutations though)
        #
        #          if there are some extra residues on ends, we can truncate the sequences
        #          maybe I wouldn't use the API, just do this alignment, if we wanted to compare structs with non-100%-identical sequences)

        alignment = next(aligner.align(seq1, seq2))
        print(alignment)
        return None

    P, Q = map(get_c_alpha_coords, (pp1, pp2))

    return rmsd.kabsch_rmsd(P, Q)


def process_group(pdb_codes):
    """ hardcoded version of the pipeline, processes a uniprot-primary group, args: pdb_codes: list

    downloads structures, filters those suitable for our analysis, subdivides by isoform, subdivides by ligand presence, runs analyses
    """

    # the following map and filter are not as readable as simple `for`s,
    #   but I think the code is this way more similar to future code, which will have configurable pipeline and operations would be parallel

    # download mmcif, build structure model, retain metadata
    structure__header__mmcif_dict = map(get_parsed_structure, pdb_codes)

    # filter structures based on model and metadata
    structures = [s for s, s_header, mmcif_dict in filter(structure_meets_our_criteria, structure__header__mmcif_dict)]

    # divide structures by isoform and ligand presence
    structs_by_isoform_by_ligand_presence = defaultdict(namedtuple('AposHolos', 'apos holos', defaults=([], [])))

    for s in structures:
        isoform = get_isoform(s)
        is_holo_ = is_holo(s, get_main_chain(s))

        isoform_group = structs_by_isoform_by_ligand_presence[isoform]

        if is_holo_:
            isoform_group.holos.append(s)
        else:
            isoform_group.apos.append(s)

    print(structs_by_isoform_by_ligand_presence)

    # analyses run for each isoform group
    for uniprot_id, isoform_group in structs_by_isoform_by_ligand_presence.items():
        # apo-holo analyses
        apo_holo_pairs = itertools.product(isoform_group.apos, isoform_group.holos)

        for apo, holo in apo_holo_pairs:
            print(apo.id, holo.id)
            print(rmsd_if_seqs_same(apo, holo))

        # holo-holo analyses
        holo_holo_pairs = itertools.combinations(isoform_group.holos, 2)

        for holo1, holo2 in holo_holo_pairs:
            print(holo1.id, holo2.id)
            print(rmsd_if_seqs_same(holo1, holo2))


def print_random_groups(n):
    """ samples uniprot groups

    Provides test input, pick suitable groups and paste them below into dict `groups`.
    Run in interactive console, so that all_uniprot_groups is loaded only once (takes long), while sampling can be run multiple times """

    all_uniprot_groups = get_basic_uniprot_groups()

    import random
    groups_sample = dict(random.sample(list(all_uniprot_groups.items()), n))
    print(groups_sample)


if __name__ == '__main__':
    # process_group(['1a1l'])

    groups = {
        # 'Q2YQQ9': defaultdict(list,
        #              {'3mqd': ['A'],
        #               '3lrf': ['A'],
        #               '4jv3': ['A'],
        #               '3u0e': ['A'],
        #               '3u0f': ['A']}),
        'A0ZZH6': defaultdict(list,
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

    for group in struct_code_groups:
        process_group(group)

    # resim pouze single-chain struktury,
    #   [kdyby byly multichain, musel bych resit 1) zda se jedna o stejny komplex (multiset uniprotů k chainům)
    #   a mapping chainů (teoreticky netriviální (víc stejných chainů), třeba superpozicí?)]
