import csv
import hashlib
import itertools
from collections import namedtuple, defaultdict
from io import StringIO
import gzip
import urllib.request
from datetime import datetime

import pandas as pd
import numpy as np
import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from pandas.core.util.hashing import hash_pandas_object, hash_array, _combine_hash_arrays


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


def get_pdb_chain_codes_with_uniprot_accession():
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz'

    with urllib.request.urlopen(url) as response:  # ~60 MB
        with gzip.open(response, 'rt', newline='', encoding='utf-8') as f:  # todo get encoding from response.headers? https://docs.python.org/3/library/urllib.request.html#urllib.response.addinfourl.headers
            date_row = next(f).strip()  # example: `# 2020/11/14 - 14:40`

            chain_mapping = pd.read_csv(f, header=0, names=(
                'pdb_code',
                'chain_id',  # in mmcif
                'uniprotkb_id',
                'label_seq_id__begin',  # in mmcif, residue sequence id
                'label_seq_id__end',
                'auth_seq_id__begin',  # in mmcif, author residue mapping to some kind of protein reference sequence, probably unimportant
                'auth_seq_id__end',
                'unp_begin',  # SIFTS mapping to the uniprotKB sequence
                'unp_end'
            ), dtype = {
                # 'pdb_code': 'category',  # s odkomentovanym vsim extra zravy a pomaly
                # 'chain_id': 'category',
                # 'uniprotkb_id': 'category',
            })

    date_csv_created = datetime.strptime(date_row.strip(' #'), '%Y/%m/%d - %H:%M')

    return chain_mapping, date_csv_created


def prepare_structure_groups(pdbchain_to_uniprot):
    gb = pdbchain_to_uniprot.groupby(['uniprotkb_id', 'unp_begin', 'unp_end'])  # or I could somehow groupby (unp_begin, unp_end) tuple if that's possible
    groups = sorted(gb, key=lambda x: len(x[1]), reverse=True)

    return groups


def _hash_dataframe_rows_no_categorize(df):
    """ returns a sequence of hashes, for each row, as in original `hash_pandas_object` """

    # hashing whole df was too slow, as it uses categorize=True for each series --> adapted pandas code

    hashes = (hash_array(series._values, categorize=False) for _, series in df.items())
    num_items = len(df.columns)
    h = _combine_hash_arrays(hashes, num_items)

    return h


def hash_dataframe_rows_no_categorize(df):
    return pd.Series(_hash_dataframe_rows_no_categorize(df), index=df.index, dtype="uint64", copy=False)


def hash_dataframe(df):
    """ returns single hash for the entire dataframe. Same hash for arbitrary row order (deliberately, by summing row-hashe) """

    return np.sum(_hash_dataframe_rows_no_categorize(df))  # deliberately summing row's hashes, allows row-order-agnostic hash


def test_string_hashing_returns_same_hash_even_if_string_are_different_in_memory():
    s1 = 'hovno'
    s2 = s1[:2] + s1[2:]

    assert s1 is not s2

    df = pd.DataFrame({'col1': [1,1], 'col2': [s1, s2]})

    row_hashes = hash_dataframe_rows_no_categorize(df)
    assert row_hashes[0] == row_hashes[1]




def equivalence_partition(iterable, relation):
    """Partitions a set of objects into equivalence classes, https://stackoverflow.com/a/38924631/3791837

    Args:
        iterable: collection of objects to be partitioned
        relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            if and only if o1 and o2 are equivalent

    Returns: classes, partitions
        classes: A sequence of sets. Each one is an equivalence class
        partitions: A dictionary mapping objects to equivalence classes
    """
    classes = []
    partitions = {}
    for o in iterable:  # for each object
        # find the class it is in
        found = False
        for c in classes:
            if relation(next(iter(c)), o):  # is it equivalent to this class?
                c.add(o)
                partitions[o] = c
                found = True
                break
        if not found:  # it is in a new class
            classes.append(set([o]))
            partitions[o] = classes[-1]
    return classes, partitions


def test():
    pdbchain_to_uniprot, _ = get_pdb_chain_codes_with_uniprot_accession()
    # struct_groups = prepare_structure_groups(pdbchain_to_uniprot)

    # chains having same segment (but how do I go to same segmentS)
    # struct_groups.filter()


    # structure ids, that have a chain (or only one) that maps to at least N residues of a uniprot id (one? because, theoretically a merged chain if there is...)
    # 1 group it by struct+chain+uniprot and count the mapping residues, filter groups of < N residues

    # pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])['uniprotkb_id', 'unp_begin', 'unp_end']


    # (determine if there's a chain with > 1 uniprot_id) Ano, napr. 6hr1
    series = pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])['uniprotkb_id'].nunique()
    series.nlargest(100)
    (series > 1).sum()  # 3584 chainů  (ale to budou asi divný struktury, beztak žádná apo+holo)

    # Now we'll inspect the uniprot segments observed in chains. For chains in an A-H group, their segment list should be the same (condition for 100 %
    # sequence identity, at least with the same length chains)
    # we'll hash the chains' dataframes, so the (rough) comparison is faster

    chain_dataframes_gb = pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])


    chain_up_map_lengths = chain_dataframes_gb.unp_end.sum() - chain_dataframes_gb.unp_begin.sum()

    # chain_dataframes = chain_dataframes_gb.filter(lambda g: chain_up_map_lengths.loc[g.name] > 50)  unfortunately, pandas filter does eager eval into a DataFrame..., losing the groups and slow

        # superslow pandas sort? When multiple columns it converts them into categoricals? Even integers?
        # .apply(lambda g: g.sort_values(['uniprotkb_id', 'unp_begin', 'unp_end'], ignore_index=True))

    #filter_condition = lambda g: chain_up_map_lengths.loc[g.name] > 50
    filter_condition = lambda g: True

    # chain_fingerprints = chain_dataframes_gb[['uniprotkb_id', 'unp_begin', 'unp_end']].apply(lambda g:
    #     hash_dataframe(g) if filter_condition(g) else None
    # ).dropna()  # unfortunate hack, pandas groupby.filter bad..,

    # made the hashing faster than above. hash rows before hashing dataframe, then .sum hashes on groupby, my hashing method already supports this
    chain_fingerprints = hash_dataframe_rows_no_categorize(pdbchain_to_uniprot[['uniprotkb_id', 'unp_begin', 'unp_end']])\
        .groupby([pdbchain_to_uniprot['pdb_code'], pdbchain_to_uniprot['chain_id']]).sum()

    chain_fingerprints = chain_fingerprints[chain_up_map_lengths > 50]  # should work even if series are in different order `pandas aligns all AXES when setting (ha- setting not indexing) Series and DataFrame from .loc, and .iloc.`

    # chain_groups_with_same_fragments = chain_fingerprints.groupby(chain_fingerprints).apply(lambda g: g.reset_index(inplace=True))

    chain_groups_with_same_fragments = [group.index.values for group_fingerprint, group in chain_fingerprints.groupby(chain_fingerprints, sort=False)]


     # does not work: sorting subdataframes too slow /sort_values(mutiple_cols) slow, that must be sorted for equals to work on dfs
     # for group_fingerprint, group in chain_fingerprints.groupby(chain_fingerprints):
     #    ideally has only one equivalence class, however if hash collision occured -- subdivide into eq. classes
     #
     #    def dataframes_equal(index_df1, index_df2):
     #        # wont work as it is not row-order agnostic (as the hash, I could sort the columns, but that I could already do in the hash--measure time)
     #        df1, df2 = map(chain_dataframes.get_group, (index_df1, index_df2))
     #        df1.reset_index(drop=True, inplace=True)  # we don't want to compare the indices (to the original df), just the contents
     #        df2.reset_index(drop=True, inplace=True)
     #        return df1.equals(df2)
     #
     #    classes, _ = equivalence_partition(group.index, dataframes_equal)
     #
     #    if len(classes) > 1:
     #        print("Lucky day")
     #
     #    chain_groups_with_same_fragments.append(classes)

    chain_groups_with_same_fragments.sort(key=lambda l: len(l), reverse=True)

    print(chain_groups_with_same_fragments[:50])
    print(len(chain_groups_with_same_fragments))

    print('groups same fragments len >1', sum(1 for _ in filter(lambda g: len(g) > 1, chain_groups_with_same_fragments)))

    struct_groups_with_multiple_structures = []

    for chg in chain_groups_with_same_fragments:
        # group by structures
        d = defaultdict(list)

        for s, ch in chg:
            d[s].append(ch)

        # aby bylo s čím porovnávat, aspon 1 apo a 1 holo => aspon 2 struktury
        if len(d) > 1:
            struct_groups_with_multiple_structures.append(d)

    print('multiple-struct-groups:', len(struct_groups_with_multiple_structures))
    print('unique_structures multiple-struct-groups', len(set(s for g in struct_groups_with_multiple_structures for s in g)))

    struct_chain_count = dict(map(lambda code: (code, 0), chain_up_map_lengths.index.get_level_values('pdb_code')))
    struct_chain_count.update(chain_up_map_lengths[chain_up_map_lengths > 50].groupby('pdb_code').count().to_dict())

    a = list(filter(
        lambda g: len(g) > 1,
        (
            list(filter(lambda struct_chains: struct_chain_count[struct_chains[0]] <= 1, g.items()))   # todo watch for == or <=
                for g in struct_groups_with_multiple_structures
        )
    ))
    print('single-chain structures: ', 'groups: ', len(a), 'unique_structures: ', len(set(s for g in a for s, chains in g)))

    print()

    # for > 50 amino uniprot mappings
    # 55444 groups same fragments len >1, 35092 actually multiple struct groups, 7124 single chain struct groups
    # into #structs:                      92375 unique structs                   36345 structures (36345 chains) (unique)

    # no minimum coverage of uniprot seq,                                       but single-chain structure has one >50 aa chain
    # 59618                               37617                                    7523
    #                                      95784                                   36770

    #                                                                           but single-chain structure has one or less >50 aa chain
    #                                                                                8285
    #                                                                               39210



    # DO MORE IN PYTHON RATHER THAN PANDAS?
    structures_with_one_chain = set()
    for pdb_code, structure_chains in pdbchain_to_uniprot.groupby(['pdb_code']):
        structure_chains.groupby('chain_id').apply(lambda r: r['unp_end'] - r['unp_begin']).filter
        # for .apply(lambda row: row.unp_end - row.unp_begin, axis=1)

test()

def test2():
    structure_code = '4cj6'
    structure_file = get_structure_stringio(structure_code)

    mmcif_parser = MMCIFParser()

    structure = mmcif_parser.get_structure(structure_code, structure_file)
    mmcif_dict_undocumented = mmcif_parser._mmcif_dict

    print(structure)
