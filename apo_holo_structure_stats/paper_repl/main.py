import itertools
import logging
import time
import warnings
from difflib import SequenceMatcher
from typing import List

import pandas as pd
from Bio.PDB import PDBList, MMCIFParser

from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure

from apo_holo_structure_stats.core.analyses import ChainResidues, ChainResidueData, ResidueId, DomainResidues
from apo_holo_structure_stats.core.analysesinstances import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.pipeline.run_analyses import chain_to_polypeptide, aligner, AnalysisHandler, JSONAnalysisSerializer

get_rotation_matrix


MIN_SUBSTRING_LENGTH_RATIO = 0.90

def compare_chains(chain1: Chain, chain2: Chain,
                   comparators__residues_param: List[Analyzer],
                   comparators__residue_ids_param: List[Analyzer],
                   comparators__domains__residues_param: List[Analyzer],
                   comparators__domains__residue_ids_param: List[Analyzer],
                   comparators__2domains__residues_param: List[Analyzer],
                   serializer_or_analysis_handler: AnalysisHandler,

) -> None:
    """ Runs comparisons between two chains. E.g. one ligand-free (apo) and another ligand-bound (holo).

    :param chain1: A Bio.PDB Chain, obtained as a part of BioPython Structure object as usual
    :param chain2: A corresponding chain (same sequence), typically from a different PDB structure. See chain1.
    """
    s1_pdb_code = chain1.get_parent().get_parent().id
    s2_pdb_code = chain2.get_parent().get_parent().id

    logging.info(f'running analyses for ({s1_pdb_code}, {s2_pdb_code}) pair...')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp1 = chain_to_polypeptide(chain1)
        pp2 = chain_to_polypeptide(chain2)

    # find longest common substring
    i1, i2, length = SequenceMatcher(a=pp1.get_sequence(), b=pp2.get_sequence(), autojunk=False).find_longest_match(0, len(pp1), 0, len(pp2))

    logging.info(f'substring to original ratio: {length/min(len(pp1), len(pp2))}')

    if length < MIN_SUBSTRING_LENGTH_RATIO * min(len(pp1), len(pp2)):
        alignment = next(aligner.align(pp1.get_sequence(), pp2.get_sequence()))
        logging.info('Sequences differ, alignment:')
        logging.info(f'\n{alignment}')
        logging.warning(f'does not meet the threshold for substring length')

    # crop polypeptides to longest common substring
    pp1 = pp1[i1:i1+length]
    pp2 = pp2[i2:i2+length]

    pp1_auth_seq_ids = {r.get_id()[1] for r in pp1}
    pp2_auth_seq_ids = {r.get_id()[1] for r in pp2}

    c1_residues = ChainResidues(list(pp1), s1_pdb_code, chain1.id)
    c2_residues = ChainResidues(list(pp2), s2_pdb_code, chain2.id)

    c1_residue_ids = ChainResidueData[ResidueId]([ResidueId.from_bio_residue(r) for r in c1_residues], s1_pdb_code, chain1.id)
    c2_residue_ids = ChainResidueData[ResidueId]([ResidueId.from_bio_residue(r) for r in c2_residues], s2_pdb_code, chain2.id)

    for a in comparators__residues_param:
        # this fn (run_analyses_for_isoform_group) does not know anything about serialization?
        # But it will know how nested it is (domain->structure) and can pass full identifiers of structures/domains

        serializer_or_analysis_handler.handle(a, a(c1_residues, c2_residues), c1_residues,
                                              c2_residues)  # in future maybe pass apo and holo. Will serialize itself. And output the object in rdf for example?
        # because what I would like is to output the analysis with objects identifiers, and then output the objects, what they contain (e.g. domain size?)


    for c in comparators__residue_ids_param:
        serializer_or_analysis_handler.handle(c, c(c1_residue_ids, c2_residue_ids), c1_residue_ids, c2_residue_ids)

    # domain-level analyses
    c1_domains = sorted(get_domains(s1_pdb_code), key=lambda d: d.domain_id)
    c2_domains = sorted(get_domains(s2_pdb_code), key=lambda d: d.domain_id)

    c1_domains__residues = [DomainResidues.from_domain(d, chain1.get_parent(), lambda id: id not in pp1_auth_seq_ids) for d in c1_domains]
    c2_domains__residues = [DomainResidues.from_domain(d, chain2.get_parent(), lambda id: id not in pp2_auth_seq_ids) for d in c2_domains]

    # assert domains equal in sequence
    # todo domain might contain the cropped residues! Fix that
    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        assert all(r1.get_resname() == r2.get_resname() for r1, r2 in zip(d_chain1, d_chain2))

    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        for a in comparators__domains__residues_param:
            serializer_or_analysis_handler.handle(a, a(d_chain1, d_chain2), d_chain1, d_chain2)

    for d_chain1, d_chain2 in zip(c1_domains, c2_domains):
        d_chain1 = d_chain1.to_set_of_residue_ids(s1_pdb_code, lambda id: id not in pp1_auth_seq_ids)
        d_chain2 = d_chain2.to_set_of_residue_ids(s2_pdb_code, lambda id: id not in pp2_auth_seq_ids)

        for a in comparators__domains__residue_ids_param:
            serializer_or_analysis_handler.handle(a, a(d_chain1, d_chain2), d_chain1, d_chain2)

    for (d1_chain1, d1_chain2), (d2_chain1, d2_chain2) in itertools.combinations(zip(c1_domains__residues, c2_domains__residues), 2):
        for a in comparators__2domains__residues_param:
            serializer_or_analysis_handler.handle(a, a(d1_chain1, d2_chain1, d1_chain2, d2_chain2), d1_chain1, d2_chain1, d1_chain2, d2_chain2)

        d1d2_chain1 = d1_chain1 + d2_chain1
        d1d2_chain2 = d1_chain2 + d2_chain2
        serializer_or_analysis_handler.handle(get_rmsd, get_rmsd(d1d2_chain1, d1d2_chain2), d1d2_chain1, d1d2_chain1)  # todo hardcoded analysis


def download_structure(pdb_code: str) -> str:
    # logging.info(f'downloading structure {pdb_code}..')

    filename = None

    for i in range(5):
        try:
            filename = PDBList(pdb='pdb_structs', verbose=False).retrieve_pdb_file(pdb_code, file_format='mmCif')
            break
        except OSError as e:
            if 'Too many connections' in str(e):
                print('too many connections')  # furt too many connections, i když stahuju přes 1 thread. Asi si to pamatuje dlouho..
                time.sleep(4)
                continue

            # structure obsolete
            try:
                logging.warning(f'{pdb_code} error: ' + str(e))
                filename = PDBList(pdb='pdb_structs', verbose=False).retrieve_pdb_file(pdb_code, file_format='mmCif', obsolete=True)  # 1seo is obsolete, but we want to replicate the paper?
                break
            except OSError as e:
                logging.warning(f'{pdb_code} error: ' + str(e))

    if filename is None:
        logging.warning(f'downloading {pdb_code} failed: ' + str(e))
        raise

    # logging.info(f'finished downloading structure {pdb_code}')
    return filename



def get_structure(pdb_code: str) -> Structure:
    filename = download_structure(pdb_code)
    return MMCIFParser(QUIET=True).get_structure(pdb_code, filename)


def get_chain_by_chain_code(s: Structure, paper_chain_code: str) -> Chain:
    if paper_chain_code == '_':
        # seems that this is the code when there's only a single chain (probably named A?)
        chain = next(iter(s[0]))

        if len(s[0]) != 1 or chain.id != 'A':
            logging.warning(f'{paper_chain_code}, underscore is not what I thought. {chain.id}, {len(s[0])}')
            return get_main_chain(s[0])  # sometime there is also chain B with ligands.
            # Get the long aa chain (probably always A, but rather use GetMainChain (with >= 50 aas))

        return chain
    return s[0][paper_chain_code]


if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)

    df = pd.read_csv('apo_holo.dat', delimiter=r'\s+', comment='#', header=None,
                     names=('apo', 'holo', 'domain_count', 'ligand_codes'), dtype={'domain_count': int})


    serializer = JSONAnalysisSerializer('output.json')

    found = False
    for index, row in df.iterrows():
        # if row.apo[:4] == '1k2x':
        #     found = True
        #
        # if not found:
        #     continue

        logging.info(f'{row.apo[:4]}, {row.holo[:4]}')

        apo = get_structure(row.apo[:4])
        holo = get_structure(row.holo[:4])

        apo_chain = get_chain_by_chain_code(apo, row.apo[4:])
        holo_chain = get_chain_by_chain_code(holo, row.holo[4:])


        compare_chains(apo_chain, holo_chain,
                       [get_rmsd, get_interdomain_surface],
                       [get_ss],
                       [get_rmsd, get_interdomain_surface],
                       [get_ss],
                       [get_hinge_angle],
                       serializer)

    serializer.dump_data()
    pass
