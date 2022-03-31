import logging
from collections import defaultdict
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from pandas._libs.parsers import STR_NA_VALUES
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from .utils import is_float, ignore_warnings

logger = logging.getLogger('scoary-load_traits')
STR_NA_VALUES.update(['-', '.'])  # compatibility with Scoary 1


class NotSplittableError(Exception):
    pass


def filter_df(df: pd.DataFrame, restrict_to: str = None, ignore: str = None) -> pd.DataFrame:
    if ignore:
        ignore = set(ignore.split(','))
        missing = ignore.difference(set(df.index))
        assert len(missing) == 0, f'Some strains in ignore were not found: {missing=}'
        df = df[[i not in ignore for i in df.index]]

    if restrict_to:
        restrict_to = set(restrict_to.split(','))
        missing = restrict_to.difference(set(df.index))
        assert len(missing) == 0, f'Some strains in restrict_to were not found: {missing=}'
        df = df[[i in restrict_to for i in df.index]]

    return df


def load_binary(traits: str, delimiter: str, restrict_to: str = None, ignore: str = None):
    dtypes = defaultdict(lambda: int)
    dtypes["index_column"] = str
    traits_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, dtype=dtypes, na_values=STR_NA_VALUES)

    filter_df(traits_df, restrict_to, ignore)
    assert traits_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert traits_df.index.is_unique, f'{traits=}: index not unique'

    # convert to bool
    traits_df = traits_df.astype('boolean')

    n_nan = int(traits_df.isna().sum().sum())
    if n_nan:
        logger.warning(f'Found {n_nan} NaN in {traits=}')

    traits_df.attrs['binarization_method'] = 'none'
    traits_df.attrs['binarization_info'] = 'none'

    return traits_df


def load_numeric(traits: str, delimiter: str):
    dtypes = defaultdict(lambda: float)
    dtypes["index_column"] = str
    numeric_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, dtype=dtypes, na_values=STR_NA_VALUES)

    assert numeric_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert numeric_df.index.is_unique, f'{traits=}: index not unique'

    n_nan = int(numeric_df.isna().sum().sum())
    if n_nan:
        logger.warning(f'Found {n_nan} NaN in {traits=}')

    return numeric_df


def apply_kmeans(kmeans: KMeans, data: pd.Series) -> (pd.Series, dict):
    nan_free_data = data[~np.isnan(data.values)]
    values = nan_free_data.values.reshape(-1, 1)
    fit = kmeans.fit(values)
    labels_ = fit.labels_.astype(bool)
    # insert_positions = [v - i for i, v in enumerate(np.where(nan_vector)[0])]
    label_to_class = dict(zip(nan_free_data.index, labels_))
    return pd.array([label_to_class.get(l, pd.NA) for l in data.index], dtype="boolean"), \
           {'method': 'kmeans', 'cutoff': np.mean(fit.cluster_centers_)}


def classify_gm(val: [float, float], certainty_cutoff: float):
    if val[0] >= certainty_cutoff:  # probability of class 1
        return True
    elif val[1] >= certainty_cutoff:  # probability of class 2
        return False
    else:
        return pd.NA


def generate_extract_covariances(covariance_type: str) -> Callable:
    def create_tuple(c1: float, c2: float) -> (float, float):
        assert type(c1) is np.float64 and type(c2) is np.float64, \
            f'Failed to read covariances! {covariance_type=} {c1=} {c2=} {type(c1)=} {type(c2)=}'
        return c1, c2

    if covariance_type == 'diag':
        return lambda covariances_: create_tuple(covariances_[0][0], covariances_[1][0])
    elif covariance_type == 'tied':
        return lambda covariances_: create_tuple(covariances_[0][0], covariances_[0][0])
    elif covariance_type == 'full':
        return lambda covariances_: create_tuple(covariances_[0][0][0], covariances_[1][0][0])
    elif covariance_type == 'spherical':
        return lambda covariances_: create_tuple(*covariances_)
    else:
        raise AssertionError('Failed to extract covariances!')


@ignore_warnings(warning=ConvergenceWarning)
def apply_gm(gm: GaussianMixture, data: pd.Series, certainty_cutoff: float, extract_covariances: Callable) -> (
        pd.Series, dict):
    nan_free_data = data[~np.isnan(data.values)]
    values = nan_free_data.values.reshape(-1, 1)
    fit = gm.fit(values)
    if not fit.converged_:
        raise NotSplittableError(f'GaussianMixture did not converge: {data.name}')
    densities = fit.predict_proba(values)
    label_to_class = {l: classify_gm(d, certainty_cutoff) for l, d in zip(nan_free_data.index, densities)}
    result = pd.array([label_to_class.get(l, pd.NA) for l in data.index], dtype="boolean")

    n_pos = result.sum()
    n_neg = (~result).sum()

    if n_pos < 1 or n_neg < 1:
        raise NotSplittableError(f'GaussianMixture failed to split {data.name}')

    metadata = {
        'method': 'gaussian', 'success': True,
        'means': [m[0] for m in gm.means_],
        'weights': list(gm.weights_),
        'covariances': extract_covariances(gm.covariances_),
    }

    return result, metadata


def parse_params(orig_params: str):
    error_message = f"""
{orig_params=} is poorly formatted.
Must be 'kmeans', 'gaussian' or '<method>:<?cutoff>:<?covariance_type>:<?alternative>'.
  Possible values for method:          {{'gaussian', 'kmeans'}}
  Possible values for delimiter:       any single character  (default: ',')
  Possible values for cutoff:          .5 <= cutoff < 1  (default: 0.85)
  Possible values for covariance_type: {{'tied', 'full', 'diag', 'spherical'}}  (default: tied)
  Possible values for alternative:     {{'skip', 'kmeans'}}  (default: skip)
""".strip()

    delimiter, cutoff, covariance_type, alternative = None, None, None, None

    params = orig_params.lower().split(':')
    method = params.pop(0)
    assert method in {'binary', 'kmeans', 'gaussian'}, error_message

    for param in params:
        if param in {'full', 'tied', 'diag', 'spherical'}:
            assert covariance_type is None, error_message
            covariance_type = param
            continue
        elif param in {'skip', 'kmeans'}:
            assert alternative is None, error_message
            alternative = param
            continue
        elif is_float(param):
            assert cutoff is None, error_message
            cutoff = float(param)
        elif len(param) == 1:
            assert delimiter is None, error_message
            delimiter = param
            continue
        else:
            raise AssertionError(error_message)

    if delimiter is None:
        delimiter = ','

    if cutoff is None:
        cutoff = 0.85

    if covariance_type is None:
        covariance_type = 'tied'

    if alternative is None:
        alternative = 'skip'

    assert 0.5 <= cutoff < 1, f'The cutoff must be between 0.5 and 1! {cutoff=}; {orig_params=}'

    return method, delimiter, cutoff, covariance_type, alternative


def binarize(
        traits: str,
        method: str, delimiter: str, cutoff: float, covariance_type: str, alternative: str,
        restrict_to: str = None,
        ignore: str = None,
        random_state: int = None,
        outfile: str = None
) -> (pd.DataFrame, pd.DataFrame):
    """
    Convert numeric data to binary data.

    :param traits: Path to traits file
    :param outfile: path to output file
    :return: numeric_df (DataFrame, dtype: float), traits_df (DataFrame, dtype: bool); columns: trait_names; index: strains
    """
    orig_numeric_df = load_numeric(traits, delimiter)

    numeric_df = filter_df(orig_numeric_df, restrict_to, ignore)

    traits_df = pd.DataFrame(index=numeric_df.index)

    kmeans = KMeans(n_clusters=2, random_state=random_state)

    binarization_info = {}

    if method == 'kmeans':
        for col in numeric_df.columns:
            binarized_data, metadata = apply_kmeans(kmeans, numeric_df[col])
            traits_df[col] = binarized_data
            binarization_info[col] = metadata

    elif method == 'gaussian':
        # apply GaussianMixture or Kmeans if GaussianMixture fails
        gm = GaussianMixture(n_components=2, random_state=random_state, covariance_type=covariance_type)
        extract_covariances = generate_extract_covariances(covariance_type)
        for col in numeric_df.columns:
            try:
                binarized_data, metadata = apply_gm(gm, numeric_df[col], certainty_cutoff=cutoff,
                                                    extract_covariances=extract_covariances)
                binarization_info[col] = metadata
            except NotSplittableError as e:
                if alternative == 'skip':
                    logger.info(f"Skipping trait '{col}': {e}.")
                    binarization_info[col] = {'method': 'gaussian', 'success': False}
                    continue
                elif alternative == 'kmeans':
                    logger.info(f"Using kmeans for trait '{col}': {e}")
                    binarized_data, metadata = apply_kmeans(kmeans, numeric_df[col])
                    binarization_info[col] = metadata
                else:
                    raise AssertionError(f'Programming error: {alternative=} must be skip or kmeans!')

            traits_df[col] = binarized_data
    else:
        raise AssertionError(f'Programming error: {method=} must be kmeans or gaussian!')

    if len(traits_df.columns) == 0:
        raise AssertionError(f'Failed to binarize traits!')

    # add metadata to DataFrame
    traits_df.attrs['binarization_method'] = method
    traits_df.attrs['binarization_info'] = binarization_info

    if outfile:
        traits_df.to_csv(outfile, sep=delimiter)

    return orig_numeric_df, traits_df


def load_traits(
        traits: str,
        trait_data_type: str = 'binary',
        restrict_to: str = None,
        ignore: str = None,
        random_state: int = None
) -> (pd.DataFrame, pd.DataFrame):
    method, delimiter, cutoff, covariance_type, alternative = parse_params(orig_params=trait_data_type)

    if method == 'binary':
        numeric_df = None
        traits_df = load_binary(
            traits=traits, delimiter=delimiter,
            restrict_to=restrict_to, ignore=ignore
        )

    else:
        numeric_df, traits_df = binarize(
            traits=traits, method=method, delimiter=delimiter, cutoff=cutoff, covariance_type=covariance_type,
            alternative=alternative, restrict_to=restrict_to, ignore=ignore, outfile=None, random_state=random_state
        )

    return numeric_df, traits_df


if __name__ == '__main__':
    import fire

    fire.Fire(binarize)
