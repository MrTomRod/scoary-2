import os
import logging
from collections import defaultdict
import json
import pandas as pd
from statsmodels.stats.multitest import multipletests
from fast_fisher.fast_fisher_numba import odds_ratio, test1t as fisher_exact_two_tailed
from queue import Empty

from .ScoaryTree import ScoaryTree
from .picking import pick
from .permutations import permute_picking
from .progressbar import print_progress
from .utils import setup_logging, AnalyzeTraitNamespace, fisher_id, grasp_namespace

logger = logging.getLogger('scoary.analyze_trait')


def worker(
        q,
        ns: AnalyzeTraitNamespace,
        result_container: {str: float | str | None},
        proc_id: int
):
    logger = setup_logging(
        logger=logging.getLogger('scoary'),
        path=f'{ns.outdir}/logs/scoary-2_proc{proc_id}.log',
        print_info=False,
        reset=True
    )
    logger.info(f'Setting up trait analysis worker {proc_id}')

    new_ns = grasp_namespace(AnalyzeTraitNamespace, ns)
    del ns

    local_result_container = {}

    while True:
        try:
            trait = q.get_nowait()
        except Empty:
            break  # completely done

        local_result_container[trait] = analyze_trait(trait, new_ns, proc_id)
        q.task_done()

    result_container.update(local_result_container)


def analyze_trait(trait: str, ns: AnalyzeTraitNamespace, proc_id: int = None) -> dict | str | None:
    logger.debug(f'Analyzing {trait=}')
    with ns.lock:
        ns.counter.value += 1
        message = trait if proc_id is None else f'P{proc_id} | {trait}'
        print_progress(
            ns.counter.value, len(ns.traits_df.columns),
            message=message, start_time=ns.start_time, message_width=25
        )

    if trait in ns.duplication_df:
        logger.debug(f'Duplicated trait: {trait} -> {ns.duplication_df[trait]}')
        save_duplicated_result(trait, ns)
        return ns.duplication_df[trait]

    trait_series = ns.traits_df[trait].dropna()
    labels = set(trait_series.index)

    if ns.all_labels == labels:
        pruned_tree = ns.tree
    else:
        pruned_tree = ns.tree.prune(labels=labels)

    result_df = init_result_df(ns.genes_bool_df, trait_series)
    test_df = create_test_df(result_df)
    test_df = add_odds_ratio(test_df)
    test_df = perform_multiple_testing_correction(
        test_df, col='fisher_p',
        method=ns.mt_f_method, cutoff=ns.mt_f_cutoff, is_sorted=True
    )

    if len(test_df) == 0:
        logger.info(f'Found 0 genes for {trait=} '
                    f'after {ns.mt_f_method}:{ns.mt_f_cutoff} filtration')
        return None

    result_df = pd.merge(
        test_df, result_df, how="left", on='__contingency_table__',
        copy=False  # for performance
    )

    result_df.attrs.update(test_df.attrs)

    if not ns.pairwise:
        result_df.attrs['best_fisher_p'] = result_df['fisher_p'].min()
        result_df.attrs['best_fisher_q'] = result_df['fisher_q'].min()
    else:
        result_df = pair_picking(
            result_df,
            significant_genes_df=ns.genes_bool_df.loc[result_df.Gene],
            tree=pruned_tree,
            label_to_trait=trait_series
        )

        if ns.worst_cutoff:
            if not (result_df['worst'] <= ns.worst_cutoff).any():
                logger.info(f'Found 0 genes for {trait=} '
                            f'after worst_cutoff={ns.worst_cutoff} filtration')
                return None

        if ns.max_genes:
            if len(result_df) > ns.max_genes:
                logger.info(f'Found too {len(result_df)} genes for {trait=} '
                            f'keeping only {ns.max_genes} with best Fisher\'s test.')
                result_df.attrs['max_genes'] = f'Trimmed {len(result_df)} genes to {ns.max_genes}.'
                result_df = result_df.iloc[:ns.max_genes]

        if ns.n_permut:
            result_df['empirical_p'] = permute_picking(
                trait=trait,
                result_df=result_df,
                tree=pruned_tree,
                label_to_trait=trait_series,
                n_permut=ns.n_permut,
                random_state=ns.random_state,
                genes_bool_df=ns.genes_bool_df
            )

            result_df['fq*ep'] = result_df['fisher_q'] * result_df['empirical_p']
            result_df.sort_values(by='fq*ep', inplace=True)

            result_df.attrs['best_fisher_p'] = result_df['fisher_p'][0]
            result_df.attrs['best_fisher_q'] = result_df['fisher_q'][0]
            result_df.attrs['best_empirical_p'] = result_df['empirical_p'][0]
            result_df.attrs['best_fq*ep'] = result_df['fq*ep'][0]

    save_result_df(trait, ns, result_df)

    # return minimal pvalues
    return result_df.attrs


def _save_trait(trait: str, ns: AnalyzeTraitNamespace):
    trait_df = pd.DataFrame(index=ns.traits_df.index)
    trait_df['class'] = ns.traits_df[trait]
    if ns.numeric_df is not None:
        trait_df['value'] = ns.numeric_df[trait]
    trait_df.index.name = 'isolate'
    trait_df.to_csv(f'{ns.outdir}/traits/{trait}/values.tsv', sep='\t')


def save_result_df(trait: str, ns: AnalyzeTraitNamespace, result_df: pd.DataFrame):
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
                 'fisher_p', 'fisher_q', 'empirical_p', 'fq*ep',
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
                info = ns.trait_info_df.loc[trait].to_dict()
                meta_data['info'] = {k: v for k, v in info.items() if not pd.isna(v)}
            except KeyError:
                pass

        json.dump(meta_data, f, indent=4, allow_nan=False)

    coverage_matrix = ns.genes_orig_df[ns.genes_orig_df.index.isin(result_df.Gene)].T
    coverage_matrix.index.name = 'Isolate'
    coverage_matrix.to_csv(f'{ns.outdir}/traits/{trait}/coverage-matrix.tsv', sep='\t')
    _save_trait(trait, ns)


def save_duplicated_result(trait: str, ns: AnalyzeTraitNamespace):
    # todo: might create symlinks traits that didn't make it through filtering!
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


def init_result_df(genes_bool_df: pd.DataFrame, trait_series: pd.Series) -> pd.DataFrame:
    """
    Create result_df with index=strains and columns=[g+t+, g+t-, g-t+, g-t-, __contingency_table__]

    :param genes_bool_df: DataFrame (dtype: bool); columns: strains; rows: genes
    :param trait_series: Boolean Series that indicates which isolates have the trait
    :return: result_df (DataFrame); columns: ['g+t+', 'g+t-', 'g-t+', 'g-t-', '__contingency_table__]; index: strains
    """
    assert trait_series.dtype == 'boolean', f'trait_series must be boolean pandas.Series!'
    assert not trait_series.hasnans, f'trait_series may not contain NANs!'
    # Preparation
    trait_pos = trait_series.index[trait_series]
    trait_neg = trait_series.index[~trait_series]
    n_pos = len(trait_pos)
    n_neg = len(trait_neg)
    n_tot = n_pos + n_neg
    assert n_tot == len(trait_series)

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
    result_df['__contingency_table__'] = [tuple(x) for x in result_df[['g+t+', 'g+t-', 'g-t+', 'g-t-']].to_numpy()]
    result_df['sensitivity'] = (result_df['g+t+'] / n_pos * 100) if n_pos else 0
    result_df['specificity'] = (result_df['g-t-'] / n_neg * 100) if n_neg else 0

    # Reset index so that Gene is its own column
    result_df.reset_index(inplace=True)

    return result_df


def create_test_df(result_df: pd.DataFrame, sort=True) -> pd.DataFrame:
    """
    Create test_df with index=__contingency_id__ and columns=[fisher_p]

    Reduce to unique contingency tables
    Add column: fisher_p

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
    test_df['fisher_p'] = test_df.__fisher_unique_table__.apply(lambda table: table_to_pval[table])

    # remove fisher_identifier
    test_df.drop('__fisher_unique_table__', axis=1, inplace=True)

    if sort:
        # sort test_df by pvalue
        test_df.sort_values(by='fisher_p', inplace=True)

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
        is_sorted=True
) -> (float, pd.DataFrame):
    assert 'fisher_p' in col
    new_col = col.replace('fisher_p', 'fisher_q')

    if method == 'native':
        reject = test_df[col] <= cutoff
        _, qval, _, _ = multipletests(pvals=test_df[col], alpha=1, method='bonferroni', is_sorted=is_sorted)
    else:
        reject, qval, alphac_sidak, alphac_bonf = multipletests(
            pvals=test_df[col],
            alpha=cutoff,
            method=method,
            is_sorted=is_sorted,
        )

    test_df[new_col] = qval
    test_df = test_df[reject]
    return test_df


def pair_picking(result_df: pd.DataFrame, significant_genes_df: pd.DataFrame, tree: ScoaryTree,
                 label_to_trait: pd.Series | dict) -> pd.DataFrame:
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
