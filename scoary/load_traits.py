import logging
import queue
from datetime import datetime
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from pandas._libs.parsers import STR_NA_VALUES
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture

from .progressbar import print_progress
from .utils import setup_logging, is_float, ignore_warnings, BinarizeTraitNamespace, MockLock, MockCounter, \
    grasp_namespace, NotSplittableError
from queue import Empty

logger = logging.getLogger('scoary.load_traits')
STR_NA_VALUES.update(['-', '.'])  # compatibility with Scoary 1


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


def switch_labels(classes: pd.Series, means: np.array) -> pd.Series:
    if means[0] > means[1]:
        return classes
    else:
        return ~classes


def load_binary(traits: str, delimiter: str, restrict_to: str = None, ignore: str = None,
                limit_traits: (int, int) = None):
    logger.debug(f'Loading binary traits: {traits=} {delimiter=}')
    traits_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, na_values=STR_NA_VALUES)

    if limit_traits is not None:
        assert len(limit_traits) == 2, f'{limit_traits=} is poorly formatted: must be (int, int).'
        l1, l2 = limit_traits
        assert type(l1) is type(l2) is int, f'{limit_traits=} is poorly formatted: must be (int, int).'
        traits_df = traits_df[traits_df.columns[l1:l2]]

    filter_df(traits_df, restrict_to, ignore)
    assert traits_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert traits_df.index.is_unique, f'{traits=}: index not unique'

    # convert to bool
    traits_df = traits_df.astype('boolean')

    n_nan = int(traits_df.isna().sum().sum())
    if n_nan:
        logger.debug(f'Found {n_nan} NaN in {traits=}')

    traits_df.attrs['binarization_method'] = 'none'
    traits_df.attrs['binarization_info'] = 'none'

    # make sure index is string and name of index is Isolate
    traits_df.index = traits_df.index.astype('str')
    traits_df.index.name = 'Isolate'

    logger.debug(f'Loaded binary traits_df:\n{traits_df}')

    return traits_df


def load_numeric(traits: str, delimiter: str, restrict_to: str = None, ignore: str = None,
                 limit_traits: (int, int) = None):
    logger.debug(f'Loading numeric traits: {traits=} {delimiter=}')
    numeric_df = pd.read_csv(traits, delimiter=delimiter, index_col=0, na_values=STR_NA_VALUES)

    if limit_traits is not None:
        assert len(limit_traits) == 2, f'{limit_traits=} is poorly formatted: must be (int, int).'
        l1, l2 = limit_traits
        assert type(l1) is type(l2) is int, f'{limit_traits=} is poorly formatted: must be (int, int).'
        numeric_df = numeric_df[numeric_df.columns[l1:l2]]

    assert numeric_df.columns.is_unique, f'{traits=}: columns are not unique'
    assert numeric_df.index.is_unique, f'{traits=}: index not unique'

    n_nan = int(numeric_df.isna().sum().sum())
    if n_nan:
        logger.debug(f'Found {n_nan} NaN in {traits=}')

    numeric_df = filter_df(numeric_df, restrict_to, ignore)

    # make sure index is string and name of index is Isolate
    numeric_df.index = numeric_df.index.astype('str')
    numeric_df.index.name = 'Isolate'

    logger.debug(f'Loaded numeric_df:\n{numeric_df}')

    return numeric_df


def apply_kmeans(kmeans: KMeans, data: pd.Series) -> (pd.Series, dict):
    nan_free_data = data[~np.isnan(data.values)]
    values = nan_free_data.values.reshape(-1, 1)
    fit = kmeans.fit(values)
    labels_ = fit.labels_.astype(bool)

    switch_labels(labels_, kmeans.cluster_centers_.flatten())

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
def apply_gm(gm: GaussianMixture, data: pd.Series, certainty_cutoff: float, extract_covariances: Callable
             ) -> (pd.Series, dict):
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

    result = switch_labels(result, gm.means_.flatten())

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
Must be 'kmeans', 'gaussian' or '<method>:<?cutoff>:<?covariance_type>:<?alternative>:<?delimiter>'.
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


def print_status(ns, trait, proc_id):
    with ns.lock:
        ns.counter.value += 1
        message = trait if proc_id is None else f'P{proc_id} | {trait}'
        print_progress(
            ns.counter.value, len(ns.numeric_df.columns),
            message=message, start_time=ns.start_time, message_width=25
        )


def fn_km(kmeans, gm, ns, trait, proc_id, container_binarized_traits, container_binarization_info, extract_covariances):
    print_status(ns, trait, proc_id)
    logger.debug(f'Binarizing {trait=} using kmeans')
    binarized_data, metadata = apply_kmeans(kmeans, ns.numeric_df[trait])
    container_binarized_traits[trait] = binarized_data
    container_binarization_info[trait] = metadata


def fn_gm(kmeans, gm, ns, trait, proc_id, container_binarized_traits, container_binarization_info, extract_covariances):
    print_status(ns, trait, proc_id)
    logger.debug(f'Binarizing {trait=} using gaussian mixture.')
    try:
        binarized_data, metadata = apply_gm(gm, ns.numeric_df[trait], certainty_cutoff=ns.cutoff,
                                            extract_covariances=extract_covariances)
        container_binarized_traits[trait] = binarized_data
        container_binarization_info[trait] = metadata

    except NotSplittableError as e:
        if ns.alternative == 'skip':
            logger.debug(f"Skipping trait '{trait}': {e}.")
            container_binarization_info[trait] = {'method': 'gaussian', 'success': False}
        elif ns.alternative == 'kmeans':
            logger.debug(f"Using kmeans for trait '{trait}': {e}")
            binarized_data, metadata = apply_kmeans(kmeans, ns.numeric_df[trait])
            container_binarized_traits[trait] = binarized_data
            container_binarization_info[trait] = metadata
        else:
            raise AssertionError(f'Programming error: {ns.alternative=} must be skip or kmeans!')


def worker(
        ns: BinarizeTraitNamespace,
        container_binarized_traits: dict = None,
        container_binarization_info: dict = None,
        queue: queue.Queue = None,
        proc_id: int = None
):
    if proc_id is None:
        # single-CPU: create vessels
        container_binarized_traits = {}
        container_binarization_info = {}
    else:
        # multiprocessing: set up logging for this process
        logger = setup_logging(
            logger=logging.getLogger('scoary'),
            path=f'{ns.outdir}/logs/scoary-2_proc{proc_id}.log',
            print_info=False,
            reset=True
        )
        logger.info(f'Setting up binarization worker {proc_id}')

    # copy all data to this Python interpreter
    new_ns = grasp_namespace(BinarizeTraitNamespace, ns)
    del ns

    kmeans = KMeans(n_init='auto', n_clusters=2, random_state=new_ns.random_state)
    gm = GaussianMixture(n_components=2, random_state=new_ns.random_state, covariance_type=new_ns.covariance_type)

    if new_ns.method == 'kmeans':
        extract_covariances = None
        fn = fn_km
    elif new_ns.method == 'gaussian':
        extract_covariances = generate_extract_covariances(new_ns.covariance_type)
        fn = fn_gm
    else:
        raise AssertionError(f'Programming error: {new_ns.method=} must be kmeans or gaussian!')

    if queue is None:
        for trait in new_ns.numeric_df.columns:
            fn(kmeans, gm, new_ns, trait, proc_id, container_binarized_traits, container_binarization_info,
               extract_covariances)
        return container_binarized_traits, container_binarization_info
    else:
        while True:
            try:
                trait = queue.get_nowait()
            except Empty:
                break  # completely done

            fn(kmeans, gm, new_ns, trait, proc_id, container_binarized_traits, container_binarization_info,
               extract_covariances)
            queue.task_done()


def binarize(
        numeric_df: pd.DataFrame,
        method: str,
        random_state: int | None,
        n_cpus: int,
        cutoff: float,
        covariance_type: str,
        alternative: str,
        outdir: str = None
) -> pd.DataFrame:
    """
    Convert numeric data to binary data.

    :param traits: Path to traits file
    :param outdir: path to outdir
    :return: numeric_df (DataFrame, dtype: float), traits_df (DataFrame, dtype: bool); columns: trait_names; index: strains
    """
    if n_cpus == 1:
        ns, lock, counter = BinarizeTraitNamespace(), MockLock(), MockCounter()
    else:
        from .init_multiprocessing import init, mp
        mgr, ns, counter, lock = init()

    logger.info(f'Binarizing traits: {method=} {alternative=} {cutoff=} {covariance_type=} {n_cpus=} {random_state=}')

    ns = BinarizeTraitNamespace.create_namespace(ns, {
        'start_time': datetime.now(),
        'lock': lock,
        'outdir': outdir,
        'counter': counter,
        'numeric_df': numeric_df,
        'method': method,
        'alternative': alternative,
        'covariance_type': covariance_type,
        'cutoff': cutoff,
        'random_state': random_state,
    })

    if n_cpus == 1:
        container_binarized_traits, container_binarization_info = worker(ns)
    else:
        mp.freeze_support()
        queue = mgr.JoinableQueue()
        container_binarized_traits = mgr.dict()
        container_binarization_info = mgr.dict()
        [queue.put(trait) for trait in numeric_df.columns]

        procs = [mp.Process(target=worker, args=(
            ns, container_binarized_traits, container_binarization_info, queue, i
        )) for i in range(n_cpus)]
        [p.start() for p in procs]
        [p.join() for p in procs]
        container_binarized_traits = dict(container_binarized_traits)
        container_binarization_info = dict(container_binarization_info)

    traits_df = pd.DataFrame.from_dict(container_binarized_traits)
    traits_df.index = numeric_df.index
    traits_df = traits_df[[c for c in numeric_df.columns if c in traits_df.columns]]

    if len(traits_df.columns) == 0:
        raise AssertionError(f'Failed to binarize traits!')

    # add metadata to DataFrame
    traits_df.attrs['binarization_method'] = method
    traits_df.attrs['binarization_info'] = container_binarization_info

    if outdir:
        traits_df.to_csv(f'{outdir}/binarized_traits.tsv', sep='\t')

    logger.debug(f'Binarized traits:\n{traits_df}')

    return traits_df


def load_traits(
        traits: str,
        trait_data_type: str = 'binary:,',
        restrict_to: str = None,
        ignore: str = None,
        n_cpus: int = 1,
        random_state: int = None,
        outdir: str = None,
        limit_traits: (int, int) = None
) -> (pd.DataFrame, pd.DataFrame):
    method, delimiter, cutoff, covariance_type, alternative = parse_params(orig_params=trait_data_type)

    if method == 'binary':
        numeric_df = None
        traits_df = load_binary(
            traits=traits, delimiter=delimiter,
            restrict_to=restrict_to, ignore=ignore,
            limit_traits=limit_traits
        )

    else:
        numeric_df = load_numeric(
            traits=traits, delimiter=delimiter,
            restrict_to=restrict_to, ignore=ignore,
            limit_traits=limit_traits
        )
        logger.info('Binarizing traits...')
        traits_df = binarize(
            numeric_df=numeric_df, method=method, random_state=random_state,
            cutoff=cutoff, covariance_type=covariance_type, alternative=alternative, n_cpus=n_cpus,
            outdir=outdir
        )
        print_progress(
            i=len(numeric_df.columns), n=len(numeric_df.columns),
            message='COMPLETE!', start_time=datetime.now(), message_width=25,
            end='\n'
        )

    return numeric_df, traits_df


if __name__ == '__main__':
    import fire

    fire.Fire(binarize)
