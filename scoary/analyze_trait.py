import os
import logging
from collections import defaultdict
import json
import pandas as pd
from statsmodels.stats.multitest import multipletests
from fast_fisher.fast_fisher_numba import odds_ratio, test1t as fisher_exact_two_tailed

from .ScoaryTree import ScoaryTree
from .picking import pick
from .permutations import permute_trait_picking
from .progressbar import print_progress
from .utils import MockNamespace, get_label_to_trait, fisher_id

logger = logging.getLogger('scoary-trait')


def analyze_trait(trait: str, ns: MockNamespace):  # -> {str: float} | str | None:
    with ns.lock:
        ns.counter.value += 1
        print_progress(
            ns.counter.value, len(ns.traits_df.columns),
            message=trait, start_time=ns.start_time, message_width=25,
            end='\n'
        )

    if trait in ns.duplication_df:
        print(f'Duplicated trait: {trait} -> {ns.duplication_df[trait]}')
        save_duplicated_result(trait, ns)
        return ns.duplication_df[trait]

    label_to_trait = get_label_to_trait(ns.traits_df[trait])

    if ns.all_labels == set(label_to_trait):
        pruned_tree = ns.tree
    else:
        pruned_tree = ns.tree.prune(labels=label_to_trait)

    result_df = init_result_df(ns.genes_bool_df, label_to_trait)
    test_df = create_test_df(result_df)
    test_df = add_odds_ratio(test_df)
    min_qval, test_df = perform_multiple_testing_correction(
        test_df, col='pval',
        method=ns.mt_f_method, cutoff=ns.mt_f_cutoff,
    )

    if len(test_df) == 0:
        logging.info(f'Found 0 genes for {trait=} '
                     f'after {ns.mt_f_method}:{ns.mt_f_cutoff} filtration')
        return None

    result_df = pd.merge(
        test_df, result_df, how="left", on='__contingency_table__',
        copy=False  # for performance
    )

    result_df.attrs.update(test_df.attrs)
    result_df.attrs['min_pval'] = result_df['pval'].min()
    result_df.attrs['min_qval'] = result_df['qval'].min()

    if not ns.no_pairwise:
        result_df = pair_picking(
            result_df,
            significant_genes_df=ns.genes_bool_df.loc[result_df.Gene],
            tree=pruned_tree,
            label_to_trait=label_to_trait
        )

        if ns.n_permut:
            result_df['pval_empirical'] = permute_trait_picking(
                result_df=result_df,
                all_label_to_gene=ns.all_label_to_gene,
                tree=pruned_tree,
                label_to_trait=label_to_trait,
                n_permut=ns.n_permut,
                random_state=ns.random_state
            )
            result_df.sort_values(by='pval_empirical', inplace=True)

            result_df.attrs['min_pval_empirical'] = result_df['pval_empirical'][0]

            # print('mtc:', ns.mt_p_method, ns.mt_p_cutoff, len(result_df))
            # result_df = perform_multiple_testing_correction(
            #     result_df, col='pval_empirical',
            #     method=ns.mt_p_method, cutoff=ns.mt_p_cutoff,
            # )
            # print('mtc:', ns.mt_p_method, ns.mt_p_cutoff, len(result_df))

            if len(result_df) == 0:
                logging.info(f'Found 0 genes for {trait=} '
                             f'after {ns.mt_p_method}:{ns.mt_p_cutoff} filtration')
                return None

            # print('permut sc1')
            # result_df['pval_empirical_scoary_1'] = calculate_confidence_interval_scoary(
            #     result_df.Gene, all_label_to_gene, tree, label_to_trait, n_permut
            # )
            # print('donepermut')

    # print(result_df.to_string())
    save_result_df(trait, ns, result_df)

    # return minimal pvalues
    return result_df.attrs


def _save_trait(trait: str, ns: MockNamespace):
    trait_df = pd.DataFrame(index=ns.traits_df.index)
    trait_df['class'] = ns.traits_df[trait]
    if ns.numeric_df is not None:
        trait_df['value'] = ns.numeric_df[trait]
    trait_df.index.name = 'isolate'
    trait_df.to_csv(f'{ns.outdir}/traits/{trait}/values.tsv', sep='\t')


def save_result_df(trait: str, ns: MockNamespace, result_df: pd.DataFrame):
    os.makedirs(f'{ns.outdir}/traits/{trait}')

    # add annotations
    if ns.gene_info_df is None:
        additional_columns = []
    else:
        additional_columns = ns.gene_info_df.columns.to_list()
        result_df = result_df.merge(ns.gene_info_df, left_on='Gene', right_index=True, how='left', copy=False)

    # reorder columns
    col_order = ['Gene'] + \
                additional_columns + \
                ['g+t+', 'g+t-', 'g-t+', 'g-t-',
                 'sensitivity', 'specificity', 'odds_ratio',
                 'pval', 'qval', 'pval_empirical',
                 'contrasting', 'supporting', 'opposing', 'best', 'worst']
    result_df = result_df[[col for col in col_order if col in result_df.columns]]

    result_df.to_csv(f'{ns.outdir}/traits/{trait}/result.tsv', sep='\t', index=False)

    binarization_info = ns.traits_df.attrs['binarization_info']
    if type(binarization_info) is str:
        binarization_info = defaultdict(lambda: 'none')

    with open(f'{ns.outdir}/traits/{trait}/meta.json', 'w') as f:
        meta_data = {
            'genes-content-type': ns.genes_orig_df.attrs['content_type'],
            'binarization-method': ns.traits_df.attrs['binarization_method'],
            'binarization-info': binarization_info[trait]
        }
        # add trait info
        if ns.trait_info_df is not None:
            try:
                meta_data['info'] = ns.trait_info_df.loc[trait].to_dict()
            except KeyError:
                meta_data['info'] = 'None'

        json.dump(meta_data, f, indent=4)

    coverage_matrix = ns.genes_orig_df[ns.genes_orig_df.index.isin(result_df.Gene)].T
    coverage_matrix.index.name = 'isolate'
    coverage_matrix.to_csv(f'{ns.outdir}/traits/{trait}/coverage-matrix.tsv', sep='\t')
    _save_trait(trait, ns)


def save_duplicated_result(trait: str, ns: MockNamespace):
    os.makedirs(f'{ns.outdir}/traits/{trait}')

    # use data from previous duplicate
    ref_trait = ns.duplication_df[trait]
    for f in ['result.tsv', 'meta.json', 'coverage-matrix.tsv']:
        os.symlink(src=f'../{ref_trait}/{f}', dst=f'{ns.outdir}/traits/{trait}/{f}')

    # create values.tsv only if numeric trait
    if ns.numeric_df is None:
        os.symlink(src=f'../{ref_trait}/values.tsv', dst=f'{ns.outdir}/traits/{trait}/values.tsv')
    else:
        _save_trait(trait, ns)


def init_result_df(genes_bool_df: pd.DataFrame, label_to_trait: dict[str:bool]) -> pd.DataFrame:
    """
    Create result_df with index=strains and columns=[g+t+, g+t-, g-t+, g-t-, __contingency_table__]

    :param genes_bool_df: DataFrame (dtype: bool); columns: strains; rows: genes
    :param trait_pos: strains that have the trait
    :param trait_neg: strains that lack the trait
    :return: result_df (DataFrame); columns: ['g+t+', 'g+t-', 'g-t+', 'g-t-', '__contingency_table__]; index: strains
    """
    # Preparation
    trait_pos = [l for l, t in label_to_trait.items() if t]
    trait_neg = [l for l, t in label_to_trait.items() if not t]
    n_pos = len(trait_pos)
    n_neg = len(trait_neg)
    n_tot = n_pos + n_neg

    # create result_df
    result_df = pd.DataFrame(index=genes_bool_df.index)
    result_df['g+t+'] = genes_bool_df[trait_pos].sum(axis=1)  # trait positive gene positive
    result_df['g+t-'] = genes_bool_df[trait_neg].sum(axis=1)  # trait negative gene positive
    result_df['g-t+'] = n_pos - result_df['g+t+']  # trait positive gene negative
    result_df['g-t-'] = n_neg - result_df['g+t-']  # trait negative gene negative

    # remove genes that are shared by none or all
    gene_sum = result_df['g+t+'] + result_df['g+t-']
    result_df = result_df[(gene_sum != 0) & (gene_sum != n_tot)]

    # Add contingency table, sensitivity and specificity
    sensitivity_fn = (lambda row: row['g+t+'] / n_pos * 100) if trait_pos else (lambda row: 0.)
    specificity_fn = (lambda row: row['g-t-'] / n_neg * 100) if trait_neg else (lambda row: 0.)

    def add_cols(row: pd.Series):
        return (
            tuple([row['g+t+'], row['g+t-'], row['g-t+'], row['g-t-']]),  # contingency table
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
    :return: test_df (DataFrame)
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


def perform_multiple_testing_correction(
        test_df: pd.DataFrame,
        method: str,
        cutoff: float,
        col: str,
) -> (float, pd.DataFrame):
    assert 'pval' in col
    new_col = col.replace('pval', 'qval')
    reject, qval, alphac_sidak, alphac_bonf = multipletests(
        pvals=test_df[col],
        alpha=cutoff,
        method=method,
        is_sorted=True,
    )
    test_df[new_col] = qval
    test_df = test_df[reject]
    min_qval = qval[0]
    return min_qval, test_df


def pair_picking(result_df: pd.DataFrame, significant_genes_df: pd.DataFrame, tree: ScoaryTree,
                 label_to_trait: {str: bool}) -> pd.DataFrame:
    """
    Required rows:
    - Gene

    Add columns:
    - Max_Pairwise_comparisons
    - Max_supporting_pairs
    - Max_opposing_pairs
    - Best_pairwise_comp_p
    - Worst_pairwise_comp_p
    """
    assert result_df.Gene.to_list() == list(significant_genes_df.index)

    max_contr, max_suppo, max_oppos, best, worst = pick(
        tree=tree.to_list, label_to_trait_a=label_to_trait,
        trait_b_df=significant_genes_df, calc_pvals=True
    )

    result_df['contrasting'] = max_contr
    result_df['supporting'] = max_suppo
    result_df['opposing'] = max_oppos
    result_df['best'] = best
    result_df['worst'] = worst

    return result_df
