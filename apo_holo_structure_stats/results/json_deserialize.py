import more_itertools
import pandas as pd


def tuple_it(l):
    if isinstance(l, list):
        return tuple(map(tuple_it, l))
    else:
        return l


def unfold_tuple_to_columns(series_or_df_with_tuple_column, new_column_names=None, column_name=None):
    """ Unfolds a column `column_name` with tuples and adds the unfolded series as new columns.
    Original column of tuples unchanged.
    Or dict.
    """
    to_unfold = series_or_df_with_tuple_column

    if isinstance(series_or_df_with_tuple_column, pd.DataFrame):
        assert column_name is not None
        to_unfold = series_or_df_with_tuple_column[column_name]

    # was too slow
    # unfolded_cols = to_unfold.apply(pd.Series)

    if isinstance(to_unfold.iloc[0], dict):
        new_column_names = new_column_names if new_column_names is not None else list(to_unfold.iloc[0].keys())
        data = list(map(list, more_itertools.unzip((tuple(d.values()) for d in to_unfold))))
    else:
        # tuple or list
        new_column_names = new_column_names if new_column_names is not None else to_unfold.name
        data = list(map(list, more_itertools.unzip(to_unfold)))

    data = dict(zip(new_column_names, data))

    return pd.DataFrame(series_or_df_with_tuple_column).assign(**data)


def tuple_columns(series1, series2):
    # the fastest way, pandas apply(tuple) reportedly 100x slower. Well,... this is Python.
    return list(zip(series1, series2))

