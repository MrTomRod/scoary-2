import pandas as pd

# from math import inf, nan, isinf, isnan, isclose

ALLOWED_CORRECTIONS = {'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by',
                       'fdr_tsbh', 'fdr_tsbky'}


def parse_correction(multiple_testing: str) -> (str, float):
    if ':' in multiple_testing:
        method, cutoff = multiple_testing.split(':', 1)
    else:
        method, cutoff = multiple_testing, 'inf'

    assert method in ALLOWED_CORRECTIONS, f'{multiple_testing=} must be in {ALLOWED_CORRECTIONS}'

    try:
        cutoff = float(cutoff)
    except ValueError:
        raise AssertionError(f'Error in {multiple_testing=}: {cutoff=} could not be converted to float')

    return method, cutoff


def get_label_to_trait(trait: pd.Series) -> {str: bool}:
    return {l: bool(t) for l, t in trait.items() if not pd.isna(t)}


def get_all_label_to_gene(genes_df: pd.DataFrame) -> pd.DataFrame:
    assert not genes_df.isna().values.any(), 'genes_df contains NaN'
    return genes_df.T.applymap(bool).to_dict()
