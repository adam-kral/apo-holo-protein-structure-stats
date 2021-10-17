import gzip
import logging
import os
import urllib.request
from collections import defaultdict, namedtuple
from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile

import requests
from requests import RequestException
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

REQUESTS_TIMEOUT = 10


class APIException(Exception):
    pass


class RequestsSessionFactory:
    retries = Retry(total=5, backoff_factor=0.3, status_forcelist=[500, 502, 503, 504])  # https://stackoverflow.com/a/35504626/3791837

    @classmethod
    def create(cls):
        s = requests.Session()
        s.mount('http://', HTTPAdapter(max_retries=cls.retries))
        s.mount('https://', HTTPAdapter(max_retries=cls.retries))
        s.headers.update({'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:87.0) Gecko/20100101 Firefox/87.0'})
        return s


def get_requests_session():
    return RequestsSessionFactory.create()


def download_and_save_file(url, file_path):
    """ Download file and return the path. Use existing file if path exists.

    Tempfile is linked to resulting file only if the download was successful, so it should be valíd. """

    r = get_requests_session().get(url, stream=True, timeout=REQUESTS_TIMEOUT)

    # if file already exists, don't download anything
    if os.path.exists(file_path):
        return  # todo won't work if file is broken, that is possible, as I added the tempfile below
        # just recently

    # download to temp file, link to actual filename only if completed
    with NamedTemporaryFile('wb', dir=os.path.dirname(file_path)) as temp:  # same dir so that os.link always works
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            temp.write(chunk)

        os.link(temp.name, file_path)  # will fail if file exists


def download_file_stringio(url):
    r = get_requests_session().get(url, stream=True, timeout=REQUESTS_TIMEOUT)

    file_like = StringIO()
    for chunk in r.iter_content(chunk_size=1024 * 1024, decode_unicode=True):
        file_like.write(chunk)

    file_like.seek(0)
    return file_like


def get_structure_stringio(code):
    return download_file_stringio(f'https://models.rcsb.org/v1/{code}/full')

def get_full_mmcif_stringio(code, dir: Path):
    url = f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{code[1:3]}/{code}.cif.gz'

    with urllib.request.urlopen(url) as response:  # ~60 MB
        with gzip.open(response, 'rt', newline='', encoding='utf-8') as r:
            with open(dir / f'{code}.cif', 'w') as f:
                f.write(r.read())


def get_best_isoform_for_chains(struct_code):
    """ Returns dict chain_id -> list of SegmentMapping(isoform_id, label_seq_id__begin, label_seq_id__end, unp_begin, unp_end, identity) """

    try:
        r = get_requests_session().get(f'https://www.ebi.ac.uk/pdbe/graph-api/mappings/isoforms/{struct_code}', timeout=REQUESTS_TIMEOUT)
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
                logging.debug(
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

    try:
        # r = get_requests_session().get(f'https://www.ebi.ac.uk/pdbe/graph-api/pdb/secondary_structure/{
        # struct_code}', timeout=REQUESTS_TIMEOUT)  # the new api endpoint was 3x slower! replaced with the old one:
        r = get_requests_session().get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/{struct_code}', timeout=REQUESTS_TIMEOUT)
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
    try:
        r = get_requests_session().get(f'https://www.ebi.ac.uk/pdbe/api/mappings/cath_b/{struct_code}', timeout=REQUESTS_TIMEOUT)
        r.raise_for_status()  # todo returns 404 if structure not found, or the data for it don't exist (e.g. tested with new release -- sifts mapping exists but there is no cath data yet)
    except RequestException as e:
        raise APIException() from e

    return r.json()[struct_code]['CATH-B']
