import random
from collections import defaultdict

from statsmodels.stats.multitest import multipletests
from scipy.stats import binom_test
import numpy as np

import logging

logger = logging.getLogger('scoary')

from .utils import *
# from .fast_fisher_numba import odds_ratio, test1t as fisher_exact_two_tailed, fisher_id
from fast_fisher.fast_fisher_cython import odds_ratio, test1t as fisher_exact_two_tailed

from .scoary_1_picking import *
import scipy.stats as ss


def fisher_id(a, b, c, d):
    """
    Eight contingency tables always give the same pvalue: ['abcd', 'acbd', 'badc', 'bdac', 'cadb', 'cdab', 'dbca', 'dcba']

    Compute and save only one version.
    """

    return min((
        (a, b, c, d),
        (a, c, b, d),
        (b, a, d, c),
        (b, d, a, c),
        (c, a, d, b),
        (c, d, a, b),
        (d, b, c, a),
        (d, c, b, a)
    ))


from .ScoaryTree import *

from pandas._libs.parsers import STR_NA_VALUES

STR_NA_VALUES.update(['-', '.'])

# from multiprocessing import Process, Manager
# _manager = Manager()
# CONFINT_CACHE = _manager.dict()
CONFINT_CACHE = {}


def scoary(
        genes: str,
        traits: str,
        start_col: int = 15,
        multiple_testing: str = 'bonferroni:0.999',
        delimiter: str = ',',
        newicktree: str = None,
        no_pairwise: bool = False,
        outdir: str = './',
        n_permut: int = 0,
        native_cutoff: float = 0.05,
        restrict_to: str = None,
        threads: int = None,
        write_tree: bool = False,
        no_time: bool = False,
        citation: bool = False,
):
    multiple_testing_method, multiple_testing_cutoff = parse_correction(multiple_testing)

    assert n_permut == 0 or n_permut >= 100, f'{n_permut=} must be at least 100.'

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

    all_labels = set(tree.labels())

    for trait in traits_df:
        print(trait)
        label_to_trait = get_label_to_trait(traits_df[trait])

        if all_labels == set(label_to_trait):
            pruned_tree = tree
        else:
            pruned_tree = tree.prune(labels=label_to_trait)

        result_df = init_result_df(genes_df, label_to_trait)
        test_df = create_test_df(result_df)
        test_df = add_odds_ratio(test_df)
        test_df = perform_multiple_testing_correction(test_df, native_cutoff, multiple_testing_method,
                                                      multiple_testing_cutoff)

        if len(test_df) == 0:
            print(f'found 0 genes for {trait=}')
            continue

        result_df = pd.merge(test_df, result_df, how="left", on='__contingency_table__',
                             copy=False)  # copy=False for performance
        # result_df = compute_pairs(result_df, all_label_to_gene, tree=pruned_tree, label_to_trait=label_to_trait)

        if n_permut:
            # conf_int = calculate_confidence_interval(genes_df, label_to_trait, n_permut=n_permut)
            # assert len(conf_int) == n_permut
            # result_df['pval_empirical'] = result_df['pval'].apply(lambda p: np.sum(conf_int <= p) / n_permut)

            result_df['pval_empirical'] = permute_acc_sc1(result_df.Gene, all_label_to_gene, tree,
                                                          label_to_trait, n_permut)

        result_df = result_df[[
            col for col in
            ('Gene', 'c1r1', 'c2r1', 'c1r2', 'c2r2', 'sensitivity', 'specificity', 'odds_ratio', 'pval', 'qval',
             'pval_empirical',
             'contrasting', 'supporting', 'opposing', 'best', 'worst'
             )
            if col in result_df.columns
        ]]
        print(result_df.to_string())


def load_genes(genes: str, delimiter: str, start_col: int) -> pd.DataFrame:
    """
    Load genes_df with columns=strains and rows=genes

    :param genes: Path to genes file
    :param delimiter: delimiter
    :param start_col: how many columns to skip
    :return: genes_df (DataFrame); columns: strains; index: genes; dtype: boolean
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

    # convert to bool
    genes_df = genes_df.astype('boolean')

    # remove core- and unique genes
    row_sums = genes_df.sum(axis=1)
    genes_df = genes_df[(row_sums != 0) & (row_sums != len(genes_df.columns))]
    return genes_df


def load_traits(traits, delimiter) -> pd.DataFrame:
    """
    Load traits_df with index=strains and columns=trait_names

    :param traits: Path to traits file
    :param delimiter: delimiter
    :return: traits_df (DataFrame); columns: trait_names; index: strains; dtype: boolean
    """
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    traits_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, dtype=dtypes, na_values=STR_NA_VALUES)
    assert traits_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert traits_df.index.is_unique, f'{traits=}: index not unique'

    # convert to bool
    traits_df = traits_df.astype('boolean')

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
    trait_pos = [l for l, t in label_to_trait.items() if t is True]
    trait_neg = [l for l, t in label_to_trait.items() if t is False]
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

    # Add contingency table, sensitivity and specificity
    sensitivity_fn = (lambda row: row['c1r1'] / n_pos * 100) if trait_pos else (lambda row: 0.)
    specificity_fn = (lambda row: row['c2r2'] / n_neg * 100) if trait_neg else (lambda row: 0.)

    def add_cols(row: pd.Series):
        return (
            tuple([row['c1r1'], row['c2r1'], row['c1r2'], row['c2r2']]),  # contingency table
            sensitivity_fn(row),
            specificity_fn(row),
        )

    result_df[['__contingency_table__', 'sensitivity', 'specificity']] = result_df.apply(
        func=add_cols, axis=1, result_type='expand')

    # Reset index so that Gene is its own column
    result_df.reset_index(inplace=True)

    return result_df


def create_test_df(result_df: pd.DataFrame, sort=True) -> pd.DataFrame:
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

    if sort:
        # sort test_df by pval
        test_df.sort_values(by='pval', inplace=True)

    return test_df


def add_odds_ratio(test_df: pd.DataFrame) -> pd.DataFrame:
    # add odds_ratio
    test_df['odds_ratio'] = test_df.__contingency_table__.apply(lambda table: odds_ratio(*table))
    return test_df


def perform_multiple_testing_correction(test_df: pd.DataFrame, native_cutoff: float, method: str,
                                        cutoff: float) -> pd.DataFrame:
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


# def compute_pairs(result_df: pd.DataFrame, all_label_to_gene, tree: ScoaryTree,
#                   label_to_trait: {str: bool}) -> pd.DataFrame:
#     """
#     Required rows:
#     - Gene
#
#     Add columns:
#     - Max_Pairwise_comparisons
#     - Max_supporting_pairs
#     - Max_opposing_pairs
#     - Best_pairwise_comp_p
#     - Worst_pairwise_comp_p
#     """
#
#     def func(row):
#         if row.odds_ratio >= 1:
#             label_to_gene = all_label_to_gene[row.Gene]
#         else:
#             label_to_gene = {l: not g for l, g in all_label_to_gene[row.Gene].items()}
#
#         contrasting = count_max_pairings(tree, label_to_trait, label_to_gene, type='contrasting')
#         # supporting = count_max_pairings(tree, label_to_trait, label_to_gene, type='supporting')
#         # opposing = count_max_pairings(tree, label_to_trait, label_to_gene, type='opposing')
#         # supporting, opposing = count_sup_op(tree, label_to_trait, label_to_gene)
#         supporting, opposing = count_best_worst(tree, label_to_trait, label_to_gene)
#
#         best = binom_test(x=supporting, n=contrasting)
#         worst = binom_test(x=contrasting - opposing, n=contrasting)
#         if worst < best:
#             best, worst = worst, best
#
#         return contrasting, supporting, opposing, best, worst
#
#     result_df[['contrasting', 'supporting', 'opposing', 'best', 'worst']] = pd.DataFrame(
#         result_df.apply(func=func, axis=1).tolist()  # autodetect dtype
#     )
#
#     return result_df


# @njit(parallel=True)
# def fisher_min_pval(test_df):
#     min_pval = 1
#
#     for i in prange(test_df.shape[0]):
#         a, b, c, d = test_df[i][0], test_df[i][1], test_df[i][2], test_df[i][3]
#
#         pval = fisher_exact_two_tailed(a, b, c, d)
#
#         if pval < min_pval:
#             min_pval = pval
#
#     return min_pval@njit(parallel=True)
def fisher_min_pval(test_df):
    min_pval = 1

    for i in range(test_df.shape[0]):
        a, b, c, d = test_df[i][0], test_df[i][1], test_df[i][2], test_df[i][3]

        pval = fisher_exact_two_tailed(a, b, c, d)

        if pval < min_pval:
            min_pval = pval

    return min_pval


def permute_acc_sc1(
        genes: pd.Series,
        all_label_to_gene: {str: bool},
        tree: ScoaryTree,
        label_to_trait: {str: bool},
        n_permut: int
) -> [int]:
    tree = tree.prune(labels=label_to_trait.keys())
    labels = tree.labels()
    tree = tree.to_list()
    boolify = lambda t1, t2: f"{'A' if t1 else 'a'}{'B' if t2 else 'b'}"
    pvals = []
    for gene in genes:
        label_to_gene = all_label_to_gene[gene]
        gtc = {l: boolify(label_to_gene[l], label_to_trait[l]) for l in labels}
        r = 0
        _MyPhyloTree_, Unpermuted_tree = convert_upgma_to_phylotree(tree, gtc)
        Pro = 'Pro' if Unpermuted_tree['Pro'] >= Unpermuted_tree['Anti'] else 'Anti'
        Unpermuted_estimator = (float(Unpermuted_tree[Pro]) / Unpermuted_tree["Total"])
        for i in range(n_permut):
            PermutedGTC = permute_gtc(gtc)
            _MyPhyloTree_, NewPhyloTree = convert_upgma_to_phylotree(tree, PermutedGTC)
            if (float(NewPhyloTree[Pro]) / NewPhyloTree["Total"] >=
                    Unpermuted_estimator):
                r += 1

            # Check how many estimators are higher than the unpermuted
            # If, after more than 30 iterations, r indicates a p > 0.1,
            # abort
            if i >= 30:
                if (1 - ss.binom.cdf(r, i, 0.1)) < 0.05:
                    emp_p = (r + 1.0) / (i + 2.0)
                    break
        else:
            emp_p = (r + 1.0) / (n_permut + 1.0)

        pvals.append(emp_p)
    return pvals


# def find_min_pval(result_df, n_pos, n_neg):
#     test_df = np.unique(result_df[['c1r1', 'c2r1', 'c1r2', 'c2r2']].to_numpy(), axis=0).astype(np.longlong)
#     min_pval = fisher_min_pval(test_df)
#     return min_pval
#
#
# def minit_result_df(genes_df: pd.DataFrame, trait_pos, trait_neg, n_tot) -> pd.DataFrame:
#     # create result_df
#     result_df = pd.DataFrame(index=genes_df.index)
#     result_df['c1r1'] = genes_df[trait_pos].sum(axis=1)  # trait positive gene positive
#     result_df['c2r1'] = genes_df[trait_neg].sum(axis=1)  # trait negative gene positive
#     result_df['c1r2'] = len(trait_pos) - result_df['c1r1']  # trait positive gene negative
#     result_df['c2r2'] = len(trait_neg) - result_df['c2r1']  # trait negative gene negative
#
#     # remove genes that are shared by none or all
#     gene_sum = result_df['c1r1'] + result_df['c2r1']
#     result_df = result_df[(gene_sum != 0) & (gene_sum != n_tot)]
#     return result_df
#
#
# def calculate_confidence_interval(genes_df: pd.DataFrame, label_to_trait: dict[str:bool], n_permut: int) -> np.array:
#     # Preparation
#     trait_pos = [l for l, t in label_to_trait.items() if t is True]
#     trait_neg = [l for l, t in label_to_trait.items() if t is False]
#     n_pos = len(trait_pos)
#     n_neg = len(trait_neg)
#
#     if (n_pos, n_neg) in CONFINT_CACHE:
#         return CONFINT_CACHE[(n_pos, n_neg)]
#
#     # filtered_genes_df = genes_df[trait_pos + trait_neg].to_numpy()
#
#     labels = list(genes_df.columns)
#     n_tot = n_pos + n_neg
#     assert n_tot <= len(labels)
#
#     with Manager() as manager:
#         shared_list = manager.list()
#
#         def permute(proc_id, shared_list) -> None:
#             permuted_labels = labels.copy()
#             random.shuffle(permuted_labels)
#
#             trait_pos = permuted_labels[:n_pos]
#             trait_neg = permuted_labels[n_neg:]
#
#             result_df = minit_result_df(genes_df, trait_pos, trait_neg, n_tot)
#             min_pval = find_min_pval(result_df, n_pos, n_neg)
#
#             shared_list.append(min_pval)
#
#         processes = [Process(target=permute, args=(i, shared_list)) for i in range(n_permut)]
#         for p in processes:
#             p.start()
#         for p in processes:
#             p.join()
#
#         result_list = np.array(shared_list)
#
#     CONFINT_CACHE[(n_pos, n_neg)] = result_list
#     return result_list


if __name__ == '__main__':
    import fire

    fire.Fire(scoary)
