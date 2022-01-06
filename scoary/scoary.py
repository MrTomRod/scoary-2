from collections import defaultdict

import pandas as pd
import numpy as np

import logging

from .tree_functions import tree_from_presence_absence, tree_from_file
from cached_contingency import CachedFisher, CachedBoschloo, odds_ratio


def scoary(
        genes: str,
        traits: str,
        start_col: int = 15,
        boschloo: bool = False,
        correction: str = 'I',
        delimiter: str = ',',
        newicktree: str = None,
        no_pairwise: bool = False,
        outdir: str = './',
        permute: bool = False,
        p_value_cutoff: float = 0.05,
        restrict_to: str = None,
        threads: int = None,
        write_tree: bool = False,
        no_time: bool = False,
        citation: bool = False,
):
    genes_df = load_genes(genes, delimiter, start_col)
    traits_df = load_traits(traits, delimiter)

    cc = CachedBoschloo() if boschloo else CachedFisher()

    # load phylogeny
    if newicktree is None:
        tree = tree_from_presence_absence(genes_df)
    else:
        with open(newicktree) as f:
            tree = tree_from_file(f.read())

    for trait in traits:
        test_df = create_test_df(traits_df.__deepcopy__(), trait)
        test_df = perform_contingency_test(test_df, cc=cc)
        test_df = compute_pairs(test_df, tree=tree)
        test_df = perform_multiple_testing_correction(test_df, correction=correction, p_value_cutoff=p_value_cutoff)


def load_genes(genes: str, delimiter: str, start_col: int) -> pd.DataFrame:
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    genes_df = pd.read_csv(genes, delimiter=delimiter, index_col=0, dtype=dtypes)
    cols_dropped = list(genes_df.columns[:start_col])
    cols_kept = list(genes_df.columns[start_col:])
    logging.info(f'Cols dropped: {cols_dropped}')
    logging.info(f'Cols kept: {cols_kept}')
    genes_df = genes_df[cols_kept]
    assert genes_df.columns.is_unique, f'{genes=}: columns not unique'
    assert genes_df.index.is_unique, f'{genes=}: index not unique'
    assert not genes_df.isna().values.any(), f'{genes=}: contains NaN'
    # remove core- and unique genes
    row_sums = genes_df.sum(axis=1)
    genes_df = genes_df[(row_sums != 0) & (row_sums != len(genes_df.columns))]
    return genes_df


def load_traits(traits, delimiter) -> pd.DataFrame:
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    traits_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, dtype=dtypes)
    assert traits_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert traits_df.index.is_unique, f'{traits=}: index not unique'
    n_nan = int(traits_df.isna().sum().sum())
    if n_nan:
        logging.warning(f'Found {n_nan} NaN in {traits=}')
    return traits_df


def create_test_df(genes_df: pd.DataFrame, trait_pos: set, trait_neg: set) -> pd.DataFrame:
    # columns: Gene, c1r1, c1r2, c2r1, c2r2, test_id, sens, spes
    n_pos = len(trait_pos)
    n_neg = len(trait_neg)
    n_tot = n_pos + n_neg

    # create test_df
    test_df = pd.DataFrame(index=genes_df.index)
    test_df['c1r1'] = genes_df[trait_pos].sum(axis=1)  # trait positive gene positive
    test_df['c2r1'] = genes_df[trait_neg].sum(axis=1)  # trait negative gene positive
    test_df['c1r2'] = n_pos - test_df['c1r1']  # trait positive gene negative
    test_df['c2r2'] = n_neg - test_df['c1r1']  # trait negative gene negative

    # remove genes that are shared by none or all
    gene_sum = test_df['c1r1'] + test_df['c2r1']
    test_df = test_df[(gene_sum != 0) & (gene_sum != n_tot)]

    print(f'number of genes: {len(test_df)}')
    return test_df


def perform_contingency_test(test_df: pd.DataFrame, cc: CachedFisher | CachedBoschloo) -> pd.DataFrame:
    # add column: pval_fisher or pval_boschloo
    # column sums are fixed
    test_df = cc.get_or_create_many(test_df=test_df)
    test_df.sort_values(by='pval', inplace=True)
    test_df.rename({'pval': f'pval_{cc.function_name}'}, axis=1, inplace=True)
    return test_df


def add_odds_ratio(test_df: pd.DataFrame) -> pd.DataFrame:
    # add column: odds_ratio
    # column sums are fixed
    test_df['odds_ratio'] = test_df.apply(lambda row: odds_ratio(row['c1r1'], row['c2r1'], row['c1r2'], row['c2r2']), axis=1)
    return test_df


def compute_pairs(data: pd.DataFrame, tree) -> pd.DataFrame:
    # add columns: Max_Pairwise_comparisons	Max_supporting_pairs	Max_opposing_pairs	Best_pairwise_comp_p	Worst_pairwise_comp_p
    return data


def perform_multiple_testing_correction(data: pd.DataFrame, correction, p_value_cutoff) -> pd.DataFrame:
    # todo: mtc based on unique contingency table or unique test string?
    return data
