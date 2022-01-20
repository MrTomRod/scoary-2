from collections import defaultdict
from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test

import logging

logger = logging.getLogger('scoary')

from fast_fisher.fast_fisher_compiled import odds_ratio, test1t as fisher_exact_two_tailed

from .utils import *
from .fisher_utils import fisher_id, contingency_id
from .ScoaryTree import ScoaryTree, count_max_pairings


def scoary(
        genes: str,
        traits: str,
        start_col: int = 15,
        multiple_testing: str = 'bonferroni:0.999',
        delimiter: str = ',',
        newicktree: str = None,
        no_pairwise: bool = False,
        outdir: str = './',
        permute: bool = False,
        native_cutoff: float = 0.05,
        restrict_to: str = None,
        threads: int = None,
        write_tree: bool = False,
        no_time: bool = False,
        citation: bool = False,
):
    multiple_testing_method, multiple_testing_cutoff = parse_correction(multiple_testing)

    # load data
    genes_df = load_genes(genes, delimiter, start_col)
    all_label_to_gene = get_all_label_to_gene(genes_df)  # {gene: {label: bool})
    traits_df = load_traits(traits, delimiter)

    # load phylogeny
    if newicktree is None:
        tree = ScoaryTree.from_presence_absence(genes_df)
    else:
        with open(newicktree) as f:
            tree = ScoaryTree.from_newick(f.read())

    for trait in traits_df:
        print(trait)
        label_to_trait = get_label_to_trait(traits_df[trait])
        result_df = init_result_df(genes_df, label_to_trait)
        test_df = create_test_df(result_df)
        test_df = add_odds_ratio(test_df)
        test_df = perform_multiple_testing_correction(test_df, native_cutoff, multiple_testing_method, multiple_testing_cutoff)
        result_df = pd.merge(test_df, result_df, how="left", on='__contingency_table__', copy=False)  # copy=False for performance
        result_df = compute_pairs(result_df, all_label_to_gene, tree=tree, label_to_trait=label_to_trait)
        print(result_df)


def load_genes(genes: str, delimiter: str, start_col: int) -> pd.DataFrame:
    """
    Load genes_df with columns=strains and rows=genes

    :param genes: Path to genes file
    :param delimiter: delimiter
    :param start_col: how many columns to skip
    :return: genes_df (DataFrame); columns: strains; index: genes; data: binary (0 or 1)
    """
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    genes_df = pd.read_csv(genes, delimiter=delimiter, index_col=0, dtype=dtypes)
    cols_dropped = list(genes_df.columns[:start_col])
    cols_kept = list(genes_df.columns[start_col:])
    logger.info(f'Cols dropped: {cols_dropped}')
    logger.info(f'Cols kept: {cols_kept}')
    genes_df = genes_df[cols_kept]
    assert genes_df.columns.is_unique, f'{genes=}: columns not unique'
    assert genes_df.index.is_unique, f'{genes=}: index not unique'
    assert not genes_df.isna().values.any(), f'{genes=}: contains NaN'
    # remove core- and unique genes
    row_sums = genes_df.sum(axis=1)
    genes_df = genes_df[(row_sums != 0) & (row_sums != len(genes_df.columns))]
    return genes_df


def load_traits(traits, delimiter) -> pd.DataFrame:
    """
    Load traits_df with index=strains and columns=trait_names

    :param traits: Path to traits file
    :param delimiter: delimiter
    :return: traits_df (DataFrame); columns: trait_names; index: strains; data: binary (0 or 1)
    """
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    traits_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, dtype=dtypes)
    assert traits_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert traits_df.index.is_unique, f'{traits=}: index not unique'
    n_nan = int(traits_df.isna().sum().sum())
    if n_nan:
        logger.warning(f'Found {n_nan} NaN in {traits=}')
    return traits_df


def init_result_df(genes_df: pd.DataFrame, label_to_trait: dict[str:bool]) -> pd.DataFrame:
    """
    Create result_df with index=strains and columns=[c1r1, c2r1, c1r2, c2r2, __contingency_table__]

    :param genes_df: DataFrame; columns: strains; rows: genes; data: binary (0 or 1)
    :param trait_pos: strains that have the trait
    :param trait_neg: strains that lack the trait
    :return: result_df (DataFrame); columns: ['c1r1', 'c2r1', 'c1r2', 'c2r2', '__contingency_table__]; index: strains
    """
    # Preparation
    trait_pos = [l for l, t in label_to_trait.items() if t]
    trait_neg = [l for l, t in label_to_trait.items() if not t]
    n_pos = len(trait_pos)
    n_neg = len(trait_neg)
    n_tot = n_pos + n_neg

    # create result_df
    result_df = pd.DataFrame(index=genes_df.index)
    result_df['c1r1'] = genes_df[trait_pos].sum(axis=1)  # trait positive gene positive
    result_df['c2r1'] = genes_df[trait_neg].sum(axis=1)  # trait negative gene positive
    result_df['c1r2'] = n_pos - result_df['c1r1']  # trait positive gene negative
    result_df['c2r2'] = n_neg - result_df['c2r1']  # trait negative gene negative

    # remove genes that are shared by none or all
    gene_sum = result_df['c1r1'] + result_df['c2r1']
    result_df = result_df[(gene_sum != 0) & (gene_sum != n_tot)]

    # add contingency table as row
    result_df['__contingency_table__'] = result_df.apply(
        func=lambda row: tuple([row['c1r1'], row['c2r1'], row['c1r2'], row['c2r2']]),
        axis=1
    )

    # reset index so that Gene is it's own column
    result_df.reset_index(inplace=True)

    return result_df


def create_test_df(result_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create test_df with index=__contingency_id__ and columns=[pval]

    Reduce to unique contingency tables
    Add column: pval_fisher

    :param result_df: DataFrame with column '__contingency_table__'
    :param sort: whether to sort the DataFrame by pvalue
    :return: test_df (DataFrame);
    """

    test_df = pd.DataFrame(result_df.__contingency_table__.unique(), columns=['__contingency_table__'])

    # add __fisher_unique_table__
    test_df['__fisher_unique_table__'] = test_df.__contingency_table__.apply(lambda table: fisher_id(*table))

    # calculate Fisher's exact test
    table_to_pval = {table: fisher_exact_two_tailed(*table) for table in test_df.__fisher_unique_table__.unique()}

    # add Fisher's exact test
    test_df['pval'] = test_df.__fisher_unique_table__.apply(lambda table: table_to_pval[table])

    # remove fisher_identifier
    test_df.drop('__fisher_unique_table__', axis=1, inplace=True)

    # sort test_df by pval
    test_df.sort_values(by='pval', inplace=True)

    return test_df


def add_odds_ratio(test_df: pd.DataFrame) -> pd.DataFrame:
    # add odds_ratio
    test_df['odds_ratio'] = test_df.__contingency_table__.apply(lambda table: odds_ratio(*table))
    return test_df


def perform_multiple_testing_correction(test_df: pd.DataFrame, native_cutoff: float, method: str, cutoff: float) -> pd.DataFrame:
    reject, qval, alphac_sidak, alphac_bonf = multipletests(
        pvals=test_df['pval'],
        alpha=cutoff,
        method=method,
        is_sorted=True,
    )
    test_df['qval'] = qval

    # filter by native pval and qval
    test_df = test_df[(test_df.pval <= native_cutoff) & reject]

    return test_df


def compute_pairs(result_df, all_label_to_gene, tree: ScoaryTree, label_to_trait: {str: bool}) -> pd.DataFrame:
    """
    Add columns:
    - Max_Pairwise_comparisons
    - Max_supporting_pairs
    - Max_opposing_pairs
    - Best_pairwise_comp_p
    - Worst_pairwise_comp_p
    """

    def func(row):
        if row.odds_ratio >= 1:
            label_to_gene = all_label_to_gene[row.Gene]
        else:
            label_to_gene = {l: not g for l, g in all_label_to_gene[row.Gene].items()}

        contrasting = count_max_pairings(tree, label_to_trait, label_to_gene, type='contrasting')
        supporting = count_max_pairings(tree, label_to_trait, label_to_gene, type='supporting')
        opposing = count_max_pairings(tree, label_to_trait, label_to_gene, type='opposing')

        best = binom_test(x=supporting, n=contrasting)
        worst = binom_test(x=contrasting - supporting, n=contrasting)
        best, worst = sorted([best, worst])

        return contrasting, supporting, opposing, best, worst

    result_df[['contrasting', 'supporting', 'opposing', 'best', 'worst']] = result_df.apply(func=func, axis=1, result_type='expand')

    return result_df


if __name__ == '__main__':
    import fire

    fire.Fire(scoary)
