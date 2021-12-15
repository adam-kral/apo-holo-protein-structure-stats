import logging
import sys
import json
from pathlib import Path
import glob

import pandas as pd

logger = logging.getLogger(__name__)


def intern_sequence_codes(items):
    """ Use string interning when parsing json AA sequence array.

    String list items (as opposed to dict keys) were not interned, thus each amino acid 3letter
    code was saved individually despite there are in total 20 or so. Took GBs and GBs of RAM
    (billions and billions,... of bytes https://youtu.be/u_aLESDql1U) when loading 300K? sequences
    """
    items['sequence'] = [sys.intern(aa_code) for aa_code in items['sequence']]
    return items


def read_json_with_seqs(path: Path) -> pd.DataFrame:
    """ Same as pd.read_json (from records), but with string interning AA 3letter codes in `sequence` column. """
    with path.open() as f:
        return pd.DataFrame.from_records(json.load(f, object_hook=intern_sequence_codes))


def maybe_print(quiet: bool, *args, **kwargs):
    if not quiet:
        print(*args, **kwargs)


def read_jsons_with_seqs_in(json_glob, quiet=False):
    """ Read multiple (or single) files specified `json_glob`.

    Show progress is quiet=False.
    """
    files = glob.glob(json_glob)
    # files = files[:10]  # for testing
    dfs = []

    for i, file in enumerate(files):
        file = Path(file)
        maybe_print(quiet, f'\rloading {file.name}: {i+1}/{len(files)}', end='')
        df = read_json_with_seqs(file)
        # drop sequence column (large)
        # some dfs are empty (low res structures in a shard)
        # if len(df):
        #     df = df.drop(columns='sequence')
        dfs.append(df)
    maybe_print(quiet)
    maybe_print(quiet, 'concatenating...')
    df = pd.concat(dfs) # concatenate all the data frames in the list.
    maybe_print(quiet, 'done.')
    return df
