import logging
from collections import defaultdict
import pandas as pd

logger = logging.getLogger('scoary.load_genes')


def filter_df(df: pd.DataFrame, restrict_to: [str] = None, ignore: [str] = None) -> pd.DataFrame:
    if ignore:
        ignore = set(ignore)
        missing = ignore.difference(set(df.columns))
        assert len(missing) == 0, f'Some strains in ignore were not found: {missing=}'
        df = df[[c for c in df.columns if c not in ignore]]

    if restrict_to is not None:
        restrict_to = set(restrict_to)
        have_cols = set(df.columns)
        cols_missing = restrict_to.difference(have_cols)
        assert len(cols_missing) == 0, f'Some strains in restrict_to were not found:' \
                                       f'\n{cols_missing=}' \
                                       f'\n{restrict_to=}' \
                                       f'\n{have_cols=}'
        cols_dropped = restrict_to.difference(set(df.columns))
        logger.debug(f'Cols kept: {list(restrict_to)}')
        logger.debug(f'Cols dropped: {list(cols_dropped)}')
        df = df[[c for c in df.columns if c in restrict_to]]

    return df


def load_gene_count_file(
        path: str,
        delimiter: str,
        restrict_to: [str] = None,
        ignore: [str] = None
) -> (pd.DataFrame, pd.DataFrame):
    """
    Load Roary-style gene count file with columns=strains and rows=genes

    :param path: Path to file
    :param delimiter: delimiter of the file
    :param restrict_to: columns to keep, will drop all other columns
    :param ignore: columns to ignore
    :return: genes_df (DataFrame, dtype: bool); columns: strains; index: genes
    """
    count_df = pd.read_csv(path, delimiter=delimiter, index_col=0)

    # remove columns that are not in traits_df
    if restrict_to is not None or ignore is not None:
        count_df = filter_df(count_df, restrict_to, ignore)

    # sanity checks
    assert count_df.columns.is_unique, f'{path=}: columns not unique'
    assert count_df.index.is_unique, f'{path=}: index not unique'
    assert not count_df.isna().values.any(), f'{path=}: contains NaN'

    # add metadata
    count_df.attrs['content_type'] = 'gene-count'

    # convert to bool
    binary_df = count_df >= 1

    # remove core- and unique genes
    row_sums = binary_df.sum(axis=1)
    binary_df = binary_df[(row_sums != 0) & (row_sums != len(binary_df.columns))]

    logger.debug(f'Loaded gene-count-df:\n{binary_df}')
    return count_df, binary_df


def load_gene_list_file(
        path: str,
        delimiter: str,
        restrict_to: [str] = None,
        ignore: [str] = None
) -> (pd.DataFrame, pd.DataFrame):
    """
    Load Orthofinder-style gene list file with columns=strains and rows=genes

    :param path: Path to file
    :param delimiter: delimiter of the file
    :param restrict_to: columns to keep, will drop all other columns
    :param ignore: columns to ignore
    :return: genes_df (DataFrame, dtype: bool); columns: strains; index: genes
    """
    list_df = pd.read_csv(path, delimiter=delimiter, index_col=0, dtype=str)

    # remove columns that are not in traits_df
    if restrict_to is not None or ignore is not None:
        list_df = filter_df(list_df, restrict_to, ignore)

    # sanity checks
    assert list_df.columns.is_unique, f'{path=}: columns not unique'
    assert list_df.index.is_unique, f'{path=}: index not unique'

    # add metadata
    list_df.attrs['content_type'] = 'gene-list'

    # convert to bool
    binary_df = ~list_df.isna()

    # remove core- and unique genes
    row_sums = binary_df.sum(axis=1)
    binary_df = binary_df[(row_sums != 0) & (row_sums != len(binary_df.columns))]

    logger.debug(f'Loaded gene-list -df:\n{binary_df}')
    return list_df, binary_df


def parse_params(orig_params: str) -> (str, str):
    error_message = f"""
{orig_params=} is poorly formatted.
Must be '<data_type>:<?delimiter>'.
  Possible values for data_type:  {{'gene-count', 'gene-list'}}  (default: gene-count)
  Possible values for delimiter:  any single character, only relevant when data_type=gene-count  (default: ',')
""".strip()

    params = orig_params.lower().split(':')

    if len(params) == 1:
        data_type, delimiter = params[0], ','
    elif len(params) == 2:
        data_type, delimiter = params
    else:
        raise AssertionError(error_message)

    assert data_type in {'gene-count', 'gene-list'}, error_message

    return data_type, delimiter


def load_genes(
        genes: str,
        gene_data_type: str,
        restrict_to: [str] = None,
        ignore: [str] = None
) -> (pd.DataFrame, pd.DataFrame):
    """
    Load genes_df with columns=strains and rows=genes

    :param genes: Path to genes file
    :return: genes_df (DataFrame, dtype: bool); columns: strains; index: genes
    """
    data_type, delimiter = parse_params(gene_data_type)

    if data_type == 'gene-count':
        genes_orig_df, genes_bool_df = load_gene_count_file(genes, delimiter, restrict_to, ignore)
    elif data_type == 'gene-list':
        genes_orig_df, genes_bool_df = load_gene_list_file(genes, delimiter, restrict_to, ignore)
    else:
        raise AssertionError(f'Programming error: {data_type=} must be gene-count or gene-list!')

    # ensure the index name is always the same
    genes_orig_df.index.name = 'Gene'
    genes_bool_df.index.name = 'Gene'

    return genes_orig_df, genes_bool_df
