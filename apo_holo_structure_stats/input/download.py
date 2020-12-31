import logging
from collections import defaultdict, namedtuple
from io import StringIO

import requests


def download_and_save_file(url, filename):
    r = requests.get(url, stream=True)

    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)


def download_file_stringio(url):
    r = requests.get(url, stream=True)

    file_like = StringIO()
    for chunk in r.iter_content(chunk_size=1024 * 1024, decode_unicode=True):
        file_like.write(chunk)

    file_like.seek(0)
    return file_like


def get_structure_stringio(code):
    return download_file_stringio(f'https://models.rcsb.org/v1/{code}/full')

# def get_full_mmcif_stringio(code):
#     url = f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{code[1:3]}/{code}.cif.gz'
#
#     with urllib.request.urlopen(url) as response:  # ~60 MB
#         with gzip.open(response, 'rt', newline='', encoding='utf-8')


def get_best_isoform_for_chains(struct_code):
    """ Returns dict chain_id -> list of SegmentMapping(isoform_id, label_seq_id__begin, label_seq_id__end, unp_begin, unp_end, identity) """

    r = requests.get(f'https://www.ebi.ac.uk/pdbe/graph-api/mappings/isoforms/{struct_code}')
    r.raise_for_status()
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

    r = requests.get(f'https://www.ebi.ac.uk/pdbe/graph-api/pdb/secondary_structure/{struct_code}')
    # if I knew entity id (not with biopython's parser), I could append it to the url

    r.raise_for_status()
    return r.json()[struct_code]['molecules']

