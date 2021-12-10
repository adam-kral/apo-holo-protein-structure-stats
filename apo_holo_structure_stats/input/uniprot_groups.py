import random
from collections import defaultdict
import gzip
import urllib.request
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
from pandas.core.util.hashing import hash_array, combine_hash_arrays

from apo_holo_structure_stats.input.download import download_and_save_file


def get_uniprot_segments_observed_dataset(download_dir=None):
    """ loads uniprot_segments_observed.csv

     the csv contains only observed segments (aa stretches with coordinates) of uniprot sequences in pdb chains"""

    # url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz'
    url = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz'

    def parse_segments_csv(file_like):
        with gzip.open(file_like, 'rt', newline='', encoding='utf-8') as f:  # todo get encoding from response.headers? https://docs.python.org/3/library/urllib.request.html#urllib.response.addinfourl.headers
            date_row = next(f).strip()  # example: `# 2020/11/14 - 14:40`

            chain_mapping_df = pd.read_csv(f, header=0, names=(
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
        return chain_mapping_df, date_csv_created

    if download_dir:
        # for reproducibility, save the file
        fname = 'uniprot_segments_observed.csv.gz'
        path = Path(download_dir) / fname
        download_and_save_file(url, path)

        with path.open('rb') as f:
            return parse_segments_csv(f)
    else:
        with urllib.request.urlopen(url) as response:  # ~60 MB
            return parse_segments_csv(response)


def get_uniprot_mappings_dataset():
    """ loads pdb_chain_uniprot.csv

    the csv contains fused observed segments (aa stretches with coordinates) together to cover the whole pdb chain sequence (incl. unobserved
    residues); the uniprot sequence maps to the pdb chain sequence even if those residues don't have coordinates observed """

    url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz'

    with urllib.request.urlopen(url) as response:  # ~60 MB
        with gzip.open(response, 'rt', newline='', encoding='utf-8') as f:  # todo get encoding from response.headers? https://docs.python.org/3/library/urllib.request.html#urllib.response.addinfourl.headers
            date_row = next(f).strip().split('|', maxsplit=1)[0]   # example: `# 2020/11/14 - 14:40 | PDB: 47.20 | UniProt: 2020.06`

            chain_mapping_df = pd.read_csv(f, header=0, names=(
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

    return chain_mapping_df, date_csv_created


def _hash_dataframe_rows_no_categorize(df):
    """ returns a sequence of hashes, for each row, as in original `hash_pandas_object` """

    # hashing whole df was too slow, as it uses categorize=True for each series --> adapted pandas code

    hashes = (hash_array(series._values, categorize=False) for _, series in df.items())
    num_items = len(df.columns)
    h = combine_hash_arrays(hashes, num_items)

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
            classes.append({o})
            partitions[o] = classes[-1]
    return classes, partitions


def get_basic_uniprot_groups(uniprot_segments_observed_df=None, download_dir=None):
    # todo neudelam to nejak rychlejsi? Problem je v tom, ze to vytvari df pro kazdou group
    # udelal bych multiindex uniprotkd_id, pdb_code, chain_id
    # velikost group muzu udelat transformem..

    # hack so we don't need more functions like this one
    if uniprot_segments_observed_df is None:
        uniprot_segments_observed_df, _ = get_uniprot_segments_observed_dataset(download_dir)

    uniprot_groupby = uniprot_segments_observed_df.drop_duplicates(['pdb_code', 'chain_id']).groupby('uniprotkb_id')[['pdb_code', 'chain_id']]
    # also possible to .groupby(['uniprotkb_id', 'unp_begin', 'unp_end']), any changes?

    uniprot_groups = {}
    for uniprotkb_id, group in uniprot_groupby:
        up_group = defaultdict(list)

        for row in group.to_records(index=False):
            up_group[row.pdb_code].append(row.chain_id)

        # we want at least 2 structures in the group (so that an apo-holo pair could exist in the group)
        if len(up_group) > 1:
            uniprot_groups[uniprotkb_id] = up_group

    return uniprot_groups


def get_basic_uniprot_groups__df(uniprot_segments_observed_df=None, download_dir=None):
    # hack so we don't need more functions like this one
    if uniprot_segments_observed_df is None:
        uniprot_segments_observed_df, _ = get_uniprot_segments_observed_dataset(download_dir)

    # todo neni jich tam vic ruznych unp segmentu k jednomu chainu (chimericky?), ty asi muzu rovnou dropnout, nebo
    # naopak nechat v obou skupinach
    df = uniprot_segments_observed_df[['uniprotkb_id', 'pdb_code', 'chain_id']]

    df = df.drop_duplicates()
    df = df.dropna()  # chain_id is NaN in extremely low number of cases, see testing_groups.ipynb (probably not, as I edited
    # this function)

    groups = df.groupby('uniprotkb_id')

    df = df.merge(groups.size().rename('uniprot_group_size'), left_on='uniprotkb_id', right_index=True)  # (need to choose one column, not important which)
    #
    # we want at least 2 structures in the group (so that an apo-holo pair could exist in the group)
    return df[df.uniprot_group_size >= 2]

def analyze_basic_uniprot_id_groups():
    """ Makes groups that consist of chains mapped to the same uniprot. This gives the upper bound of
    structures/chains/groups (and can estimate than the number of pairs) we would analyse
    """
    uniprot_segments_observed_df, _ = get_uniprot_segments_observed_dataset()

    chains = get_basic_uniprot_groups__df(uniprot_segments_observed_df)
    # chains = chains.set_index(['uniprotkb_id', 'pdb_code', 'chain_id'])

    print('multiple-struct-groups:', chains['uniprotkb_id'].nunique())
    print('unique_structures multiple-struct-groups', chains['pdb_code'].nunique())
    print(chains.iloc[:5])

    print(chains.describe())

    # chain_dataframes_gb = uniprot_segments_observed_df.groupby(['pdb_code', 'chain_id'])
    # chain_up_map_lengths = chain_dataframes_gb.unp_end.sum() - chain_dataframes_gb.unp_begin.sum()
    #
    # struct_chain_count = dict(map(lambda code: (code, 0), chain_up_map_lengths.index.get_level_values('pdb_code')))
    # struct_chain_count.update(chain_up_map_lengths[chain_up_map_lengths > 15].groupby('pdb_code').count().to_dict())
    #
    # single_chain_structs = list(filter(
    #     lambda g: len(g) > 1,
    #     (
    #         list(filter(lambda struct_chains: struct_chain_count[struct_chains[0]] <= 1, g.items()))  # todo watch for == or <=
    #         for g in uniprot_groups.values()
    #     )
    # ))
    #
    # print('single-chain structures:: ',
    #       'groups: ', len(single_chain_structs),
    #       'unique_structures: ', len(set(s for g in single_chain_structs for s, chains in g)))
    # multiple-struct-groups: 27092
        # unique_structures multiple-struct-groups 139913
    # single-chain (one chain > 50) structures:  groups:  9159 unique_structures:  60294
    # single-chain (one chain > 15) structures:  groups:  8683 unique_structures:  57308

    # #chains in a group, usually < 30, rarely > 100



# použít místo uniprot_segments_observed.csv.gz tohle: pdb_chain_uniprot.csv.gz, tam jsou extendlý ty alignmenty, a pokud je na jeden chain víc fragmentů od jednoho uniprotu,
# znamená to rozdíl v sekvenci (ale i jen v SEQRES, nemusí být ten rozdíl observed).
def analyze_uniprot_groups_joined_fragments():
    """ Makes groups that require all uniprot segments mapped to a chain to be equal, for each chain in a group.
    That probably was not the case of the structures in the paper. Structures had leading/trailing residues, along
    the LCS which was used for the analysis. These trailing residues would probably be in uniprot, the segments
    wouldn't be equal, so they wouldn't be in the group. So this estimate is lower than 'analyze_basic_uniprot_id_groups',
    but still the actual number of groups/pairs will probably be lower (due to filtering criteria like resolution, ligand-free
    and ligand-bound in one group,...)
    """
    pdbchain_to_uniprot, _ = get_uniprot_segments_observed_dataset()

    # pdbchain_to_uniprot, _ = get_uniprot_mappings_dataset()1111

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
    print('multiple-struct-groups total len', sum(len(g) for g in struct_groups_with_multiple_structures))

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

    ## uniprot_segments_observed.csv

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


    ## pdb_chain_uniprot.csv
    # groups same fragments len >1 54328
    # multiple-struct-groups: 31126  # proč je těch prvních dvou míň? Protože tam možná můžou být mutace a moje groupovani segmentu neni tak prisne (nezna sekvence, jen indexy), jako je algoritmus sifts
    # unique_structures multiple-struct-groups 122018
    # single-chain structures:  groups:  8214 unique_structures:  47936


if __name__ == '__main__':
    analyze_basic_uniprot_id_groups()
    # analyze_uniprot_groups_joined_fragments()
