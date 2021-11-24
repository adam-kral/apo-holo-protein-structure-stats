import pandas as pd


def tuple_it(l):
    if isinstance(l, list):
        return tuple(map(tuple_it, l))
    else:
        return l


def unfold_tuple_to_columns(series_or_df_with_tuple_column, new_column_names=None, column_name=None):
    """ Unfolds a column `column_name` with tuples and adds the unfolded series as new columns.
    Original column of tuples unchanged. """
    to_unfold = series_or_df_with_tuple_column

    if isinstance(series_or_df_with_tuple_column, pd.DataFrame):
        assert column_name is not None
        to_unfold = series_or_df_with_tuple_column[column_name]

    unfolded_cols = to_unfold.apply(pd.Series)  # todo why this works? What if the cols don't have the same dtype?
    if new_column_names:
        unfolded_cols.columns = new_column_names
    return pd.concat([series_or_df_with_tuple_column, unfolded_cols], axis=1)


def tuple_columns(series1, series2):
    # the fastest way, pandas apply(tuple) reportedly 100x slower. Well,... this is Python.
    return list(zip(series1, series2))


