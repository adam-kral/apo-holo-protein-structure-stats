from collections import defaultdict
import gzip
import urllib.request
from datetime import datetime

import pandas as pd
import numpy as np
from pandas.core.util.hashing import hash_array, _combine_hash_arrays


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
                # 'pdb_code': 'category',  # when anything uncommented, extra memory and time hungry, category probably suited only for a small set of values
                # 'chain_id': 'category',
                # 'uniprotkb_id': 'category',
            })

    date_csv_created = datetime.strptime(date_row.strip(' #'), '%Y/%m/%d - %H:%M')

    return chain_mapping, date_csv_created


def _hash_dataframe_rows_no_categorize(df):
    """ returns a sequence of hashes, for each row, as in original `hash_pandas_object` """

    # hashing whole df was too slow, as it uses categorize=True for each series --> adapted pandas code

    hashes = (hash_array(series._values, categorize=False) for _, series in df.items())
    num_items = len(df.columns)
    h = _combine_hash_arrays(hashes, num_items)

    return h


def hash_dataframe_rows_no_categorize(df):
    """ returns a Series of hashes, for each row, as in original `hash_pandas_object` """

    return pd.Series(_hash_dataframe_rows_no_categorize(df), index=df.index, dtype="uint64", copy=False)


def hash_dataframe(df):
    """ returns single hash for the entire dataframe. Same hash for any row ordering (deliberately, by _summing_ row-hashes) """

    return np.sum(_hash_dataframe_rows_no_categorize(df))  # deliberately summing row's hashes, allows row-order-agnostic hash


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


def analyze_basic_uniprot_id_groups():
    pdbchain_to_uniprot, _ = get_pdb_chain_codes_with_uniprot_accession()
    uniprot_groupby = pdbchain_to_uniprot.drop_duplicates(['pdb_code', 'chain_id']).groupby('uniprotkb_id')[['pdb_code', 'chain_id']]
    # also possible to .groupby(['uniprotkb_id', 'unp_begin', 'unp_end'])

    uniprot_groups = {}
    for uniprotkb_id, group in uniprot_groupby:
        up_group = defaultdict(list)

        for row in group.to_records(index=False):
            up_group[row.pdb_code].append(row.chain_id)

        if len(up_group) > 1:
            uniprot_groups[uniprotkb_id] = up_group

    print('multiple-struct-groups:', len(uniprot_groups))
    print('unique_structures multiple-struct-groups', len(set(s for g in uniprot_groups.values() for s in g)))

    chain_dataframes_gb = pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])
    chain_up_map_lengths = chain_dataframes_gb.unp_end.sum() - chain_dataframes_gb.unp_begin.sum()

    struct_chain_count = dict(map(lambda code: (code, 0), chain_up_map_lengths.index.get_level_values('pdb_code')))
    struct_chain_count.update(chain_up_map_lengths[chain_up_map_lengths > 15].groupby('pdb_code').count().to_dict())

    single_chain_structs = list(filter(
        lambda g: len(g) > 1,
        (
            list(filter(lambda struct_chains: struct_chain_count[struct_chains[0]] <= 1, g.items()))  # todo watch for == or <=
            for g in uniprot_groups.values()
        )
    ))

    print('single-chain structures:: ',
          'groups: ', len(single_chain_structs),
          'unique_structures: ', len(set(s for g in single_chain_structs for s, chains in g)))
    # multiple-struct-groups: 27092
    # unique_structures multiple-struct-groups 139913
    # single-chain (one chain > 50) structures:  groups:  9159 unique_structures:  60294
    # single-chain (one chain > 15) structures:  groups:  8683 unique_structures:  57308

analyze_basic_uniprot_id_groups()


# tohle je možná trochu blbost (ale sežrala čas), observed fragmenty u různých struktur stejných proteinů mohou být různé -- v jednom experimentu něco být observed nemusí, trochu jinej krystal
# ALE todo použít místo uniprot_segments_observed.csv.gz tohle: pdb_chain_uniprot.csv.gz, tam jsou extendlý ty alignmenty, a pokud je na jeden chain víc fragmentů od jednoho uniprotu,
# znamená to rozdíl v sekvenci (ale i jen v SEQRES, nemusí být ten rozdíl observed).
def analyze_uniprot_groups_all_fragments_identical():
    pdbchain_to_uniprot, _ = get_pdb_chain_codes_with_uniprot_accession()

    # chains having same segment (but how do I go to same segmentS)


    # structure ids, that have a chain (or only one) that maps to at least N residues of a uniprot id (one? because, theoretically a merged chain if there is...)
    # 1 group it by struct+chain+uniprot and count the mapping residues, filter groups of < N residues



    # (determine if there's a chain with > 1 uniprot_id) Ano, napr. 6hr1, v dokumentaci `chimeric`
    series = pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])['uniprotkb_id'].nunique()
    print((series > 1).sum())  # 3584 chainů  (ale to budou asi divný struktury, beztak žádná apo+holo)

    # Now we'll inspect the uniprot segments observed in chains (structure). For chains in an A-H group, their segment list should (hopefully, see input
    # docs) be the same (condition for 100 % sequence identity, at least with the same length chains)

    # we'll hash the chains' dataframes, so the (rough) comparison is faster

    chain_dataframes_gb = pdbchain_to_uniprot.groupby(['pdb_code', 'chain_id'])

    # made the hashing faster than before. hash rows before hashing dataframe, then .sum hashes on groupby, my hashing concept already supports this
    chain_fingerprints = hash_dataframe_rows_no_categorize(pdbchain_to_uniprot[['uniprotkb_id', 'unp_begin', 'unp_end']])\
        .groupby([pdbchain_to_uniprot['pdb_code'], pdbchain_to_uniprot['chain_id']]).sum()

    chain_up_map_lengths = chain_dataframes_gb.unp_end.sum() - chain_dataframes_gb.unp_begin.sum()

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
     #        print("Lucky day, hash collision")
     #
     #    chain_groups_with_same_fragments.append(classes)

    # have the largest groups of chains first
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
    struct_chain_count.update(chain_up_map_lengths[chain_up_map_lengths > 15].groupby('pdb_code').count().to_dict())  # now 15

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

    #                                                                            single-chain structure has one or less >15-aa chain
    #                                                                             6934
    #                                                                             35299

    # no minimum coverage of uniprot seq...                                    ...but single-chain structure has one >50-aa chain (allows e.g.  40-aa chain appearing in the groups, larger than 1st)
    # 59618                               37617                                    7523
    #                                      95784                                   36770
    #                                                                           having <= 1 >50-aa chain
    #                                                                           7652
    #                                                                           37002

analyze_uniprot_groups_all_fragments_identical()
