import os
import json
import logging
from collections import defaultdict
import numpy as np
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
        step: int,
        result_container: {dict | str | None},
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

    analyze_trait_fn = analyze_trait_step_1_fisher if step == 1 else analyze_trait_step_2_pairpicking

    local_result_container = {}

    while True:
        try:
            trait = q.get_nowait()
        except Empty:
            break  # completely done

        local_result_container[trait] = analyze_trait_fn(trait, new_ns, proc_id)
        q.task_done()

    result_container.update(local_result_container)


def analyze_trait_step_1_fisher(trait: str, ns: AnalyzeTraitNamespace, proc_id: int = None) -> np.ndarray | str:
    logger.debug(f"Analyzing {trait=}, step 1: Fisher's test")
    with ns.lock:
        ns.counter.value += 1
        message = trait if proc_id is None else f'P{proc_id} | {trait}'
        print_progress(
            ns.counter.value, ns.queue_size,
            message=message, start_time=ns.start_time, message_width=25
        )

    if trait in ns.duplicates:
        logger.debug(f'Duplicated trait: {trait} -> {ns.duplicates[trait]}')
        save_duplicated_result(trait, ns)
        return ns.duplicates[trait]

    # Prepare results.tsv
    isolate_trait_series = ns.traits_df[trait].dropna()
    result_df = init_result_df(ns.genes_bool_df, isolate_trait_series)

    # Sometimes, binarization gives extreme results and no genes are left
    if len(result_df) == 0:
        logger.info(f'Found 0 genes for {trait=}!')
        return False

    # Compute Fisher's test efficiently
    test_df = create_test_df(result_df)
    test_df = add_odds_ratio(test_df)
    result_df = pd.merge(test_df, result_df, how="left", on='__contingency_table__', copy=False)

    # Perform multiple testing correction
    multiple_testing_df = result_df[['__pattern_id__', 'fisher_p']].drop_duplicates('__pattern_id__')
    if ns.trait_wise_correction:
        multiple_testing_df = multiple_testing_correction(
            multiple_testing_df, 'fisher_p', 'fisher_q',
            ns.mt_f_method, ns.mt_f_cutoff, True
        )
        if len(multiple_testing_df) == 0:
            logger.info(f'Found 0 genes for {trait=} after multiple testing correction!')
            return False

        multiple_testing_df.drop('fisher_p', axis=1, inplace=True)
        result_df = pd.merge(multiple_testing_df, result_df, how="left", on='__pattern_id__', copy=False)
        result = True
    else:
        result = multiple_testing_df

    os.makedirs(f'{ns.outdir}/traits/{trait}')
    result_df.to_csv(f'{ns.outdir}/traits/{trait}/result.tsv', sep='\t', index=False)

    return result


def analyze_trait_step_2_pairpicking(trait: str, ns: AnalyzeTraitNamespace, proc_id: int = None) -> dict | str | None:
    logger.debug(f'Analyzing {trait=}, step 2: Pair picking')
    with ns.lock:
        ns.counter.value += 1
        message = trait if proc_id is None else f'P{proc_id} | {trait}'
        print_progress(
            ns.counter.value, ns.queue_size,
            message=message, start_time=ns.start_time, message_width=25
        )
    summary_data = {}

    result_df = pd.read_csv(f'{ns.outdir}/traits/{trait}/result.tsv', sep='\t')

    if ns.trait_wise_correction:
        assert 'fisher_q' in result_df.columns, f'{result_df.columns=} must contain "fisher_q"!'
    else:
        multiple_testing_df = ns.multiple_testing_df.loc[trait, :]
        result_df = pd.merge(multiple_testing_df, result_df, how="left", on='__pattern_id__', copy=False)

    assert 'fisher_p' in result_df.columns, f'{result_df.columns=} must contain "fisher_p"!'
    assert 'fisher_q' in result_df.columns, f'{result_df.columns=} must contain "fisher_q"!'

    if not ns.pairwise:
        min_row = result_df.loc[result_df['fisher_p'].idxmin()]
        summary_data['best_fisher_p'] = min_row['fisher_p']
        summary_data['best_fisher_q'] = min_row['fisher_q']
    else:
        trait_series = ns.traits_df[trait].dropna()
        isolates = set(trait_series.index)
        if ns.all_labels == isolates:
            pruned_tree = ns.tree
        else:
            pruned_tree = ns.tree.prune(labels=isolates)

        result_df = pair_picking(
            result_df,
            significant_genes_df=ns.genes_bool_df.loc[result_df.Gene],
            tree=pruned_tree,
            label_to_trait=trait_series
        )

        if ns.worst_cutoff:
            keep = result_df['worst'] <= ns.worst_cutoff
            if not keep.any():
                logger.info(f'Found 0 genes for {trait=} '
                            f'after worst_cutoff={ns.worst_cutoff} filtration')
                return None
            result_df = result_df[keep]

        assert result_df.fisher_p.is_monotonic_increasing, f'{result_df.fisher_p=} must be monotonic increasing!'

        if ns.max_genes:
            if len(result_df) > ns.max_genes:
                logger.info(f'Found too {len(result_df)} genes for {trait=} '
                            f'keeping only {ns.max_genes} with best Fisher\'s test.')
                summary_data['max_genes'] = f'Trimmed {len(result_df)} genes to {ns.max_genes}.'
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

            best_row = result_df.iloc[0]
            summary_data['best_fisher_p'] = best_row['fisher_p']
            summary_data['best_fisher_q'] = best_row['fisher_q']
            summary_data['best_empirical_p'] = best_row['empirical_p']
            summary_data['best_fq*ep'] = best_row['fq*ep']

    save_result_df(trait, ns, result_df)

    # return minimal pvalues
    return summary_data


def _save_trait(trait: str, ns: AnalyzeTraitNamespace):
    trait_df = pd.DataFrame(index=ns.traits_df.index)
    trait_df['binary'] = ns.traits_df[trait]
    if ns.numeric_df is not None:
        trait_df['numeric'] = ns.numeric_df[trait]
    trait_df.index.name = 'isolate'
    trait_df.to_csv(f'{ns.outdir}/traits/{trait}/values.tsv', sep='\t')


def save_result_df(trait: str, ns: AnalyzeTraitNamespace, result_df: pd.DataFrame):
    # add annotations
    if ns.gene_info_df is None:
        additional_columns = []
    else:
        additional_columns = ns.gene_info_df.columns.to_list()
        result_df = result_df.merge(ns.gene_info_df, left_on='Gene', right_index=True, how='left', copy=False)

    # reorder columns
    col_order = ['Gene', *additional_columns,
                 'g+t+', 'g+t-', 'g-t+', 'g-t-',
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
    os.makedirs(f'{ns.outdir}/traits/{trait}')

    # use data from previous duplicate
    ref_trait = ns.duplicates[trait]
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

    # Create result_df
    result_df = pd.DataFrame(index=genes_bool_df.index)
    result_df['g+t+'] = genes_bool_df[trait_pos].sum(axis=1)  # trait positive gene positive
    result_df['g+t-'] = genes_bool_df[trait_neg].sum(axis=1)  # trait negative gene positive
    result_df['g-t+'] = n_pos - result_df['g+t+']  # trait positive gene negative
    result_df['g-t-'] = n_neg - result_df['g+t-']  # trait negative gene negative

    # Remove genes that are shared by none or all
    gene_sum = result_df['g+t+'] + result_df['g+t-']
    to_keep = (gene_sum != 0) & (gene_sum != n_tot)
    result_df = result_df[to_keep]

    # Add unique pattern ID
    genes_bool_df_reduced = genes_bool_df.loc[to_keep, trait_pos.to_list() + trait_neg.to_list()]
    pattern_id = genes_bool_df_reduced.groupby(by=genes_bool_df_reduced.columns.to_list()).ngroup()
    result_df['__pattern_id__'] = pattern_id

    # Add contingency table, sensitivity and specificity
    result_df['__contingency_table__'] = [tuple(x) for x in result_df[['g+t+', 'g+t-', 'g-t+', 'g-t-']].to_numpy()]
    if n_pos:
        pos_sensitivity = (result_df['g+t+'] / n_pos * 100)  # use if positive g/t correlation
        neg_sensitivity = (result_df['g-t+'] / n_pos * 100)  # use if negative g/t correlation
    else:
        pos_sensitivity = neg_sensitivity = pd.Series(0, index=result_df.index)

    if n_neg:
        pos_specificity = (result_df['g-t-'] / n_neg * 100)  # use if positive g/t correlation
        neg_specificity = (result_df['g+t-'] / n_neg * 100)  # use if negative g/t correlation
    else:
        pos_specificity = neg_specificity = pd.Series(0, index=result_df.index)

    keep_pos = (pos_sensitivity + pos_specificity) > (neg_sensitivity + neg_specificity)
    result_df["sensitivity"] = pos_sensitivity.where(keep_pos, neg_sensitivity)
    result_df["specificity"] = pos_specificity.where(keep_pos, neg_specificity)

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


def multiple_testing_correction(
        df: pd.DataFrame,
        pval_column: str,
        qval_column: str,
        method: str,
        cutoff: float,
        is_sorted: bool = False
) -> (float, pd.DataFrame):
    assert pval_column in df.columns, f'{pval_column=} must be in {df.columns=}!'
    if qval_column in df.columns:
        logger.warning(f'Overwriting {qval_column=} in {df.columns=}!')

    pvals = df[pval_column]

    # Apply multiple testing correction for each orthogene
    if method == 'native':
        reject = pvals <= cutoff
        _, qval, _, _ = multipletests(pvals=pvals, alpha=1, method='bonferroni', is_sorted=is_sorted)
    else:
        reject, qval, alphac_sidak, alphac_bonf = multipletests(
            pvals=pvals, alpha=cutoff, method=method, is_sorted=is_sorted,
        )

    df[qval_column] = qval
    df = df[reject]
    return df


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
