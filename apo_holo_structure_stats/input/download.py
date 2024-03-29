import gzip
import logging
import os
import urllib.request
import warnings
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Set, Union, Tuple, overload, Literal, Dict, Any

import requests
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Structure import Structure
from requests import RequestException
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

from apo_holo_structure_stats.settings import Settings
from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds

logger = logging.getLogger(__name__)

REQUESTS_TIMEOUT = Settings.API_REQUESTS_RETRIES


class APIException(Exception):
    pass


class RequestsSessionFactory:
    """ This Session retries requests on errors. """
    # https://stackoverflow.com/a/35504626/3791837
    retries = Retry(total=Settings.API_REQUESTS_RETRIES, backoff_factor=0.3, status_forcelist=[500, 502, 503, 504])

    @classmethod
    def create(cls):
        s = requests.Session()
        s.mount('http://', HTTPAdapter(max_retries=cls.retries))
        s.mount('https://', HTTPAdapter(max_retries=cls.retries))
        # todo why :)?
        s.headers.update({'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:87.0) Gecko/20100101 Firefox/87.0'})
        return s


def get_requests_session():
    return RequestsSessionFactory.create()


def download_and_save_file(url, file_path):
    """ Download file and return the path.

    Tempfile is linked to resulting file only if the download was successful, so it should be valíd.

    :param url:
    :param file_path:
    """

    # download to temp file, link to actual filename only if completed
    with NamedTemporaryFile('wb', dir=os.path.dirname(file_path)) as temp:  # same dir so that os.link always works

        with get_requests_session() as s:
            for chunk in s.get(url, stream=True, timeout=REQUESTS_TIMEOUT).iter_content(chunk_size=1024 * 1024):
                temp.write(chunk)

        if os.path.exists(file_path):
            os.remove(file_path)

        os.link(temp.name, file_path)  # will fail if file exists
        # os.replace wouldn't work easily, as the exit context manager for tempfile will crash (No such file or dir..)

#
# def download_file_stringio(url):
#     r = get_requests_session().get(url, stream=True, timeout=REQUESTS_TIMEOUT)
#
#     file_like = StringIO()
#     for chunk in r.iter_content(chunk_size=1024 * 1024, decode_unicode=True):
#         file_like.write(chunk)
#
#     file_like.seek(0)
#     return file_like
#
#
# def get_structure_stringio(code):
#     return download_file_stringio(f'https://models.rcsb.org/v1/{code}/full')
#
# def get_full_mmcif_stringio(code, dir: Path):
#     url = f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{code[1:3]}/{code}.cif.gz'
#
#     with urllib.request.urlopen(url) as response:  # ~60 MB
#         with gzip.open(response, 'rt', newline='', encoding='utf-8') as r:
#             with open(dir / f'{code}.cif', 'w') as f:
#                 f.write(r.read())


def get_best_isoform_for_chains(struct_code):
    """ Get best isoform from pdbe API

    (however, the results were not correct AFAIK - sometimes the isoform had worse
    aligment than the perfect alignment with another isoform (the actual best was the primary accession, but reported as
    best was some isoform accession)

    Returns dict chain_id -> list of SegmentMapping(isoform_id, label_seq_id__begin, label_seq_id__end, unp_begin, unp_end, identity) """

    with get_requests_session() as s:
        try:
            r = s.get(f'https://www.ebi.ac.uk/pdbe/graph-api/mappings/isoforms/{struct_code}', timeout=REQUESTS_TIMEOUT)
            r.raise_for_status()
        except RequestException as e:
            raise APIException() from e

    data = r.json()[struct_code]['UniProt']

    chain_based_data = defaultdict(list)

    SegmentMapping = namedtuple('SegmentMapping', 'uniprot_id label_seq_id__begin label_seq_id__end unp_begin unp_end identity')

    for isoform_id, isoform_data in data.items():
        for mapping in isoform_data['mappings']:
            chain_based_data[mapping['chain_id']].append(  # chain_id is probably? auth_asym_id and  struct_asym_id is I would say label_asym_id
                SegmentMapping(
                    isoform_id,
                    mapping['start']['residue_number'],  # teď nevim co by bylo mapping['pdb_start'], viz debug print nize
                    mapping['end']['residue_number'],
                    mapping['unp_start'],
                    mapping['unp_end'],
                    mapping['identity'],
                ))

            # exploratory debug print, for development
            if (mapping['start']['residue_number'], mapping['end']['residue_number']) != (mapping['pdb_start'], mapping['pdb_end']):
                logger.debug(
                    f'pdbe API response, Get Best Isoform for {struct_code}:  `pdb_` sequence positions do not match those in `start` and `end`')

    return chain_based_data


def get_secondary_structure(struct_code: str):
    """ Get secondary structure info from pdbe API

    (that is based on their computation and is consistent across entries unlike SS info in mmcif files)

    :return: list of molecule objects (loaded json at "molecules")

    docs: https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/#api-PDB-GetPDBSecondaryStructures

    Example API json response:
    {
        "1cbs": {
            "molecules": [
                {
                    "entity_id": 1,
                    "chains": [
                        {
                            "chain_id": "A",
                            "struct_asym_id": "A",
                            "secondary_structure": {
                                "helices": [
                                    {
                                        "start": {
                                            "author_residue_number": 14,
                                            "author_insertion_code": null,
                                            "residue_number": 14
                                        },
                                        "end": {
                                            "author_residue_number": 22,
                                            "author_insertion_code": null,
                                            "residue_number": 22
                                        }
                                    }
                                ],
                                "strands": [
                                    {
                                        "start": {
                                            "author_residue_number": 5,
                                            "author_insertion_code": null,
                                            "residue_number": 5
                                        },
                                        "end": {
                                            "author_residue_number": 13,
                                            "author_insertion_code": null,
                                            "residue_number": 13
                                        },
                                        "sheet_id": 1
                                    }
                                ]
                            }
                        }
                    ]
                }
            ]
        }
    }
    """
    with get_requests_session() as s:
        try:
            # r = s.get(f'https://www.ebi.ac.uk/pdbe/graph-api/pdb/secondary_structure/{
            # struct_code}', timeout=REQUESTS_TIMEOUT)  # the new api endpoint was 3x slower! replaced with the old one:
            r = s.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/{struct_code}', timeout=REQUESTS_TIMEOUT)
            # if I knew entity id (not with biopython's parser), I could append it to the url
            r.raise_for_status()
        except RequestException as e:
            raise APIException() from e

    return r.json()[struct_code]['molecules']


def get_domains(struct_code: str):
    """ Get CATH-B domains from pdbe API.

    docs: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html#sifts_apidiv_call_5_calltitle

    How are the domains represented?

    The domains are composed of residue segments: CATH's API has pdb_segments array: http://www.cathdb.info/version/v4_1_0/api/rest/domain_summary/1cukA01
    Domains are not always contiguous. A domain can be inserted in between residues belonging to another one (e.g. 1ex6 in tests). The
    inserted residues (domain) then form their own structurally packed domain.

    More on CATH domain representation: http://wiki.cathdb.info/wiki-beta/doku.php?id=data:cathdomall

    See the example response below and the return statement for the return value of this function.

    Example API json response:
    {
      "1cbs": {
        "CATH-B": {
          "2.40.128.20": {
            "homology": "Calycin beta-barrel core domain",
            "mappings": [
              {
                "domain": "1cbsA00",
                "end": {
                  "author_residue_number": 137,
                  "author_insertion_code": "",
                  "residue_number": 137
                },
                "segment_id": 1,
                "entity_id": 1,
                "chain_id": "A",
                "start": {
                  "author_residue_number": 1,
                  "author_insertion_code": "",
                  "residue_number": 1
                },
                "struct_asym_id": "A"
              }
            ],
            "topology": "Lipocalin",
            "architecture": "Beta Barrel",
            "identifier": "Lipocalin",
            "class": "Mainly Beta",
            "name": "Cellular retinoic acid binding protein type ii. Chain: a. Engineered:yes"
          }
        }
      }
    }

    Could not use newer API https://www.ebi.ac.uk/pdbe/graph-api/mappings/cath/1bs4, as for this particular structure for example
    it returns wrong domain end residue numbers - 168 instead of 1168 -  so it reports the segment is 1001 - 168. As of 7/23/2021.
    TODO!! Už to snad spravili, můžu použít nový. Ne je to nějak divně, místo 168, 500něco, 1168 tam je furt 1168
    """
    with get_requests_session() as s:
        try:
            r = s.get(f'https://www.ebi.ac.uk/pdbe/api/mappings/cath_b/{struct_code}', timeout=REQUESTS_TIMEOUT)
            r.raise_for_status()  # todo returns 404 if structure not found, or the data for it don't exist (e.g. tested
            #   with new release -- sifts mapping exists but there is no cath data yet)
        except RequestException as e:
            raise APIException() from e

    return r.json()[struct_code]['CATH-B']


def find_or_download_structure(pdb_code: str, allow_download=True) -> Path:
    """ Download structure file and return its path. Use existing file, if path already exists. """
    filename = f'{pdb_code}.cif.gz'
    url = f'https://files.rcsb.org/download/{filename}'

    local_path = Settings.STRUCTURE_STORAGE_DIRECTORY / filename

    # if file already exists, don't download anything
    # else download it
    if os.path.exists(local_path):
        logger.info(f'using cached structure {pdb_code} at {local_path}')
    elif allow_download:
        # create dir structure if does not exist
        Settings.STRUCTURE_STORAGE_DIRECTORY.mkdir(parents=True, exist_ok=True)

        download_and_save_file(url, local_path)
        logger.info(f'finished downloading structure {pdb_code} to {local_path}')
    else:
        raise ValueError('Structure not found on disk and download not allowed')

    return local_path


class CustomMMCIFParser(MMCIFParser):
    """ Adapted BioPython code, just to get the pdb code (_entry.id) from the mmcif file """
    def get_structure(self, structure_id, file=None, mmcif_dict=None):
        """ Parses file contents and returns Structure object.

        Note that parameter order is different to the BioPython's implementation (reversed, as structure_id is optional).

        :param structure_id: if None, taken from mmcif (`_entry.id`)
        :param file: a file-like object or a file name
        :param mmcif_dict: this may be specified instead of file to build the structure from already parsed mmcif dict
        :return: Bio.PDB.Structure
        """

        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            # begin change
            self._mmcif_dict = MMCIF2Dict(file) if not mmcif_dict else mmcif_dict

            if structure_id is None:
                structure_id = self._mmcif_dict['_entry.id'][0].lower()
            # end change

            self._build_structure(structure_id)
            self._structure_builder.set_header(self._get_header())

        return self._structure_builder.get_structure()


@dataclass
class MmcifParseResult:
    structure: Structure # Biopython's Structure object
    bio_to_mmcif_mappings: BiopythonToMmcifResidueIds.Models  # maps auth_seq_ids (used in BioPython) to
    # label_seq_ids, a modern way to identify a residue
    poly_seqs: BiopythonToMmcifResidueIds.EntityPolySequences  # polymer sequences, as in mmcif
    poly_with_microhet: Set[int]  # ids of sequences with micro-heterogeneity (multiple resolved AAs for a position)


# maybetodo refactor this, allow just parsing the mmcif and then getting the structure subsequently (so that the SMCRA model
# does not have to be built for structures that don't meet basic params (e.g. resolution). I hope and suspect that
# building the model would take the majority of time compared to just parsing the mmcif into a dict extremely large
# 3J3Q took 2 min of cpu time and approx 3 gb mem to load the mmcif dict and mmcifparser.get_structure took ..5 GB
# memory and 3.5 cpu time, so ... don't do it, I only reserved too few memory (however with the caching -- this might
# be a problem if there will be some large structures next to each other in 4 consecutive structures (cache size 3).
# And there are --> 3j34, 3j3q, 3j3y so... that won't fit the memory unless I have ~16 gb per node

 # bipythoní parser neni ideální, využívá legacy PDB fieldy (ATOM/HETATM) a auth_seq_id (s fallbackem na label_seq_id == problém při mapování z pdbe api), auth_asym_id. Pokud by např.
    # nebyl u HETATM auth_seq_id (nepovinný ale, in about 100.0 % of entries), spadlo by to
    # auth_seq_id u heteroatomů (v mmcifu mají všechny heteroatomy label_seq_id `.`) umožnuje identifikaci, do jaké molekuly atom patří (jinak by byl jen název sloučeniny)
    # ovšem ty auth_ položky nemusí být číslem, ale např. tento parser je převádí na int() -- může spadnout
    #    zde https://bioinformatics.stackexchange.com/a/11590 jsem se ale dočetl undocumented fakt,
    #    že ty písmenka nějak dali pryč a převedli je do insertion code sloupce
    # Anyway, I need label_seq_id of Residue objects not only for correct mapping of pdbe API results (see above), but also to align whole
    # polypeptide sequence (even unobserved residues) to form valid chain (apo-holo) pairs.



def parse_mmcif(pdb_code: str = None, path: Union[str, os.PathLike] = None, with_mmcif_dict=False,
                allow_download=True) -> Union[MmcifParseResult, Tuple[MmcifParseResult, MMCIF2Dict]]:
    """ Parses mmcif file and builds the biopython SMCRA model (structure, model, chain, residue, atom) of the structure.

    pdb_code can be specified, in that case Settings.STRUCTURE_STORAGE_DIRECTORY is inspected (where ah-download-structures
    writes).
    If path is specified, structures is loaded from that file. File can be gzipped (.gz) or plain text (arbitrary or no
    extension).
    :param pdb_code:
    :param path:
    :param with_mmcif_dict: whole mmcif dict, suspect may be large, therefore this class is only returned when needed
    :param allow_download:
    :return:
    """
    local_path = find_or_download_structure(pdb_code, allow_download) if not path else path

    text_file = gzip.open(local_path, 'rt', newline='', encoding='utf-8')\
        if str(local_path).endswith('.gz') else open(local_path, 'rt', newline='', encoding='utf-8')

    with text_file:
        mmcif_parser = CustomMMCIFParser(QUIET=True)  # todo quiet good idea?

        structure = mmcif_parser.get_structure(pdb_code, text_file)
        # reuse already parsed mmcifdict (albeit undocumented)
        mapping, poly_seqs, with_microhet = BiopythonToMmcifResidueIds.create(mmcif_parser._mmcif_dict)

        result = MmcifParseResult(structure, mapping, poly_seqs, with_microhet)

        if with_mmcif_dict:
            return result, mmcif_parser._mmcif_dict
        return result
