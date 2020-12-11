import itertools
from collections import defaultdict, namedtuple
import logging
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser, is_aa, NeighborSearch, PDBList


from apo_holo_structure_stats.input.download import get_structure_stringio, get_best_isoform_for_chains
from apo_holo_structure_stats.input.uniprot_groups import get_basic_uniprot_groups


def get_parsed_structure(structure_code):
    structure_file = PDBList(pdb=pathlib.Path.cwd() / 'pdb_structs').retrieve_pdb_file(structure_code, file_format='mmCif')  # (mmCif is the default, but implicit = warning)

    # bipythoní parser neni ideální, využívá legacy PDB fieldy (ATOM/HETATM) a auth_seq_id (s fallbackem na label_seq_id == problém při mapování z pdbe api), auth_asym_id. Pokud by např.
    # nebyl u HETATM auth_seq_id (nepovinný ale, in about 100.0 % of entries), spadlo by to
    #  auth_asym_id, in about 93.8 % of entries -- zde to může spadnout
    # auth_seq_id u heteroatomů (v mmcifu mají všechny heteroatomy label_seq_id `.`) umožnuje identifikaci, do jaké molekuly atom patří (jinak by byl jen název sloučeniny)
    # ovšem ty auth_ položky nemusí být číslem, ale např. tento parser je převádí na int() -- může spadnout
    mmcif_parser = MMCIFParser()

    structure = mmcif_parser.get_structure(structure_code, structure_file)
    mmcif_dict_undocumented = mmcif_parser._mmcif_dict

    return structure, mmcif_parser.header, mmcif_dict_undocumented


def single_struct_analyzer(struct):
    pass


def get_hetero_atom_residues(struct):
    """ non-polymer ligands, excl. water

    residue BioPython's meaning is a group of atoms belonging to same
    scans all chains (we deliberately ignore author chain assignment in case of a heteroatom)"""

    return filter(lambda res: res.id[0].startswith('H_'), struct.get_residues())  # .id[0] is a residue's hetatm flag


def get_short_peptide_ligands(struct, peptide_length_limit):
    return filter(lambda chain: sum(is_aa(residue) for residue in chain) <= peptide_length_limit, struct.get_chains())


def get_ligands(struct):
    return itertools.chain(get_hetero_atom_residues(struct), get_short_peptide_ligands(struct, 15))


def get_chains(struct):
    return filter(lambda chain: sum(is_aa(residue) for residue in chain) >= 50, struct.get_chains())


def is_holo(struct, main_chain):
    def has_at_least_n_non_hydrogen_atoms(ligand, n):
        non_hydrogen_atoms = 0

        for atom in ligand.get_atoms():
            assert atom.element is not None
            if atom.element != 'H':
                non_hydrogen_atoms += 1

            if non_hydrogen_atoms >= n:
                return True   # todo nakonec můžu asi sumovat všechny, stejně budu chtít konfigurovatelny output, aby mi dal počet atomů ligandu, nebo budu dělat statistiky, kolik atomu má průměrný ligand atp.


        return False

    # ligand has >= 6 non-hydrogen atoms
    ligands = list(filter(lambda lig: has_at_least_n_non_hydrogen_atoms(lig, 6), get_ligands(struct)))

    # ligand is within RADIUS in contact with MIN_RESIDUES_WITHIN_LIGAND residues

    # (in original paper they used a program LPC, for ensuring specific interaction of ligand with at least 6 residue, this is a "shortcut",
    #    a temporary condition (simple))

    main_chain_atoms = list(main_chain.get_atoms())
    ns = NeighborSearch(main_chain_atoms)

    RADIUS = 4.5
    MIN_RESIDUES_WITHIN_LIGAND = 4

    acceptable_ligands = []

    for ligand in ligands:
        residues_in_contact_with_ligand = set()
        for ligand_atom in ligand.get_atoms():  # ligand can be a chain or a residue

            main_chain_atoms_in_contact = ns.search(ligand_atom.get_coord(), RADIUS)

            for atom in main_chain_atoms_in_contact:
                residues_in_contact_with_ligand.add(atom.get_parent())

        if len(residues_in_contact_with_ligand) >= MIN_RESIDUES_WITHIN_LIGAND:
            acceptable_ligands.append(ligand)

    return len(acceptable_ligands) > 0


def structure_meets_our_criteria(s, s_header, mmcif_dict):
    """ decides if structure meets criteria for resolution, single-chainedness, etc. """

    resolution = s_header['resolution']
    print(resolution)

    # skip low resolution
    if resolution and resolution > 2.5:
        return False  # todo skip structure, asi aby se nastavil reason: filtered out: resolution > 2.5

    # skip non-xray. Done in the original paper. Try without it. Or then, specify in output dataset which experimental
    #   method (or one could integrate the info themselves)
    pass

    # skip DNA/RNA complexes
    try:
        if any(poly_type in mmcif_dict['_entity_poly.type'] for poly_type in ('polydeoxyribonucleotide', 'polyribonucleotide')):
            return False
    except KeyError:
        # _entity_poly.type:: Required in PDB entries no; Used in current PDB entries Yes, in about 100.0 % of entries
        # could also do the check manually (i.e. exists a chain with >= 2 nucleotides, but I wouldn't know for sure if they form a polymer)
        #       in the paper they allow only single nucleotides
        pass

    # skip non-single-chain structure (too short chain or too many chains)
    chains_for_analysis = list(get_chains(s))
    if len(chains_for_analysis) != 1:
        return False

    return True


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
        #          if there are some extra residues on ends, we can truncate the sequences (if not, skip, or in future, truncate)
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

    provides test input, pick suitable groups and paste them below into dict `groups`

    run in interactive console, so that all_uniprot_groups is loaded only once (takes long), while sampling can be run multiple times """

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