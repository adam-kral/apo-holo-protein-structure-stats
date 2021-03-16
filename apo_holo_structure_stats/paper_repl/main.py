import itertools
import logging
import time
import warnings
from difflib import SequenceMatcher
from typing import List

import pandas as pd
from Bio.PDB import PDBList, MMCIFParser

from Bio.PDB.Chain import Chain
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.Structure import Structure

from apo_holo_structure_stats.core.analyses import ChainResidues, ChainResidueData, ResidueId, DomainResidues, DomainResidueData
from apo_holo_structure_stats.core.analysesinstances import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.input.download import APIException
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

    # get domains (set of auth_seq_id), sort them by domain id and hope they will correspond to each other
    # or could map corresponding domains by choosing the ones that have the most overlap?
    try:
        c1_domains = sorted(filter(lambda d: d.chain_id == chain1.id, get_domains(s1_pdb_code)), key=lambda d: d.domain_id)
        c2_domains = sorted(filter(lambda d: d.chain_id == chain2.id, get_domains(s2_pdb_code)), key=lambda d: d.domain_id)
    except APIException as e:
        if e.__cause__ and '404' in str(e.__cause__):
            logging.warning('no domains found, skip the domain-level analysis')
            return  # no domains found, skip the domain-level analysis
        raise


    assert len(c1_domains) == len(c2_domains)

    # SequenceMatcher on domain resiudes
    c1_domains__residues = []
    c2_domains__residues = []

    for c1_d, c2_d in zip(c1_domains, c2_domains):
        c1_d_residues = DomainResidues.from_domain(c1_d, chain1.get_parent(), lambda id: id not in pp1_auth_seq_ids)
        c2_d_residues = DomainResidues.from_domain(c2_d, chain2.get_parent(), lambda id: id not in pp2_auth_seq_ids)

        c1_d_pp = Polypeptide(c1_d_residues)
        c2_d_pp = Polypeptide(c2_d_residues)

        # find longest common substring of domains
        i1, i2, length = SequenceMatcher(a=c1_d_pp.get_sequence(), b=c2_d_pp.get_sequence(), autojunk=False).find_longest_match(0, len(c1_d_pp), 0,
                                                                                                                        len(c2_d_pp))
        logging.info(f'domain substring to original ratio: {length / min(len(c1_d_pp), len(c2_d_pp))}')

        if length < MIN_SUBSTRING_LENGTH_RATIO * min(len(c1_d_pp), len(c2_d_pp)):
            alignment = next(aligner.align(c1_d_pp.get_sequence(), c2_d_pp.get_sequence()))
            logging.info('Sequences differ, alignment:')
            logging.info(f'\n{alignment}')
            logging.warning(f'does not meet the threshold for substring length')

        # todo reflect domain cropping in the object id (domain id) somehow?
        c1_domains__residues.append(DomainResidues(c1_d_residues[i1:i1+length], c1_d_residues.structure_id, c1_d_residues.chain_id, c1_d_residues.domain_id))
        c2_domains__residues.append(DomainResidues(c2_d_residues[i2:i2+length], c2_d_residues.structure_id, c2_d_residues.chain_id, c2_d_residues.domain_id))

    # assert domains equal in sequence
    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        assert len(d_chain1) == len(d_chain2)
        assert all(r1.get_resname() == r2.get_resname() for r1, r2 in zip(d_chain1, d_chain2))

    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        for a in comparators__domains__residues_param:
            serializer_or_analysis_handler.handle(a, a(d_chain1, d_chain2), d_chain1, d_chain2)

    for d_chain1, d_chain2 in zip(c1_domains__residues, c2_domains__residues):
        # Convert DomainResidues to DomainResidueData[ResidueId]
        d_chain1 = DomainResidueData[ResidueId]([ResidueId.from_bio_residue(r) for r in d_chain1], d_chain1.structure_id, d_chain1.chain_id, d_chain1.domain_id)
        d_chain2 = DomainResidueData[ResidueId]([ResidueId.from_bio_residue(r) for r in d_chain2], d_chain2.structure_id, d_chain2.chain_id, d_chain2.domain_id)

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
        if row.apo[:4] == '1cgj':
            continue  # chymotrypsinogen x chymotrypsin + obojí má ligand... (to 'apo' má 53 aa inhibitor)

        if row.apo[:4] == '1ehd':
            found = True

        if not found:
            continue

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
