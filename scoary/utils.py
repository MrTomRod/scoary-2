import importlib.metadata
import os
import sys
import json
import logging
from copy import deepcopy
import warnings
from functools import cache
from typing import Type, Any, Callable
from datetime import datetime
import numpy as np
import pandas as pd
from importlib.metadata import version

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
ALLOWED_CORRECTIONS = {'native', 'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh',
                       'fdr_by', 'fdr_tsbh', 'fdr_tsbky'}

logger = logging.getLogger('scoary.utils')

try:
    from ete3 import Tree as EteTree


    def print_tree(scoary_tree, label_to_gene: {str: bool}, label_to_trait: {str: bool}, show_label=True):
        if show_label:
            label_fn = lambda label: f'{int(label_to_gene[label])}{int(label_to_trait[label])}_{label}'
        else:
            label_fn = lambda label: f'{int(label_to_gene[label])}{int(label_to_trait[label])}'

        renamed_tree = scoary_tree.rename(label_fn)
        ete_tree = EteTree(renamed_tree.to_newick())
        print(ete_tree)

except ImportError as e:
    def print_tree(scoary_tree, label_to_gene: {str: bool}, label_to_trait: {str: bool}):
        raise ImportError('This function requires the ete3 library. Please install via "pip install ete3"')


@cache
def get_version() -> str:
    try:
        return version('scoary-2')
    except importlib.metadata.PackageNotFoundError:
        return 'development'


class NotSplittableError(Exception):
    pass


class NoTraitsLeftException(Exception):
    pass


def decode_unicode(string: str) -> str:
    return string.encode('utf-8').decode('unicode-escape')


def setup_outdir(outdir: str, input: dict) -> str:
    outdir = outdir.rstrip('/')
    assert not os.path.exists(outdir), f'ERROR: {outdir=} already exists!'
    os.makedirs(f'{outdir}/traits')
    os.makedirs(f'{outdir}/logs')
    os.makedirs(f'{outdir}/app')
    with open(f'{outdir}/logs/input.json', 'w') as f:
        json.dump(input, f, indent=4)
    return outdir


def setup_logging(logger: logging.Logger, path: str = None, print_info: bool = True, reset: bool = False):
    """
    Setup logging for Scoary

    :param logger: logging.logging.Logger
    :param path: if set, DEBUG and higher goes to log files
    :param print_info: if True, INFO and higher goes to stdout
    :param reset: if True: close and remove all file handlers. (Important for multiprocessing: removes locks!)
    :return:
    """
    if reset:
        while logger.handlers:
            handler = logger.handlers[0]
            handler.close()
            logger.removeHandler(handler)

    logger.setLevel(logging.DEBUG)

    if path is not None:
        # create logfile
        logfile = logging.FileHandler(path)
        logfile.setLevel(logging.DEBUG)
        logfile.setFormatter(logging.Formatter("%(asctime)s [%(name)s: %(levelname)s] %(message)s"))
        logger.addHandler(logfile)

    if print_info:
        # create streamhandler
        stdout = logging.StreamHandler()
        stdout.setLevel(logging.INFO)
        logger.addHandler(stdout)

    return logger


def ignore_warnings(warning: Type[Warning]):
    """
    Decorator to suppress warnings.

    Example:

    @ignore_warnings(warning=ConvergenceWarning)
    def some_function():
        # any produced ConvergenceWarnings will be suppressed

    :param warning: class of warning to be suppressed
    """

    def decorator(function):
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', warning)
                return function(*args, **kwargs)

        return wrapper

    return decorator


class RecursionLimit:
    """
    Context manager to temporarily change the recursion limit.

    Example:

    with RecursionLimit(10 ** 6):
        some_function_that_requires_lots_of_recursions()
    """
    new: int
    old: int

    def __init__(self, new_recursion_limit: int):
        self.new = new_recursion_limit

    def __enter__(self):
        self.old = sys.getrecursionlimit()
        logger.debug(f'Setting new recursion limit: {self.old} -> {self.new}')
        sys.setrecursionlimit(self.new)

    def __exit__(self, *args, **kwargs):
        logger.debug(f'Setting old recursion limit: {self.new} -> {self.old}')
        sys.setrecursionlimit(self.old)


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


def is_int(string: str) -> bool:
    try:
        int(string)
        return True
    except ValueError:
        return False


def is_float(string: str) -> bool:
    try:
        float(string)
        return True
    except ValueError:
        return False


def split_into_parts(list_: list, n_parts: int) -> [list]:
    quotient, reminder = divmod(len(list_), n_parts)
    return [
        list_[i * quotient + min(i, reminder):(i + 1) * quotient + min(i + 1, reminder)]
        for i in range(n_parts)
    ]


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


def load_info_file(
        logger: logging.Logger,
        info_file: str,
        merge_col: str,
        expected_overlap_set: set = None,
        reference_file: str = None
) -> pd.DataFrame:
    """
    Load an info_file into pd.DataFrame:
        - Separator: tab ('\t')
        - Must have this header: {merge_col}\t{colname1}\t{colname2}...

    :param logger: instance of logging.Logger
    :param info_file: path to file
    :param merge_col: name of first column
    :param expected_overlap_set: a set of strings, some of which must occur in the index of info_file
    :param reference_file: path to reference file, just used for error messages
    :return: pd.DataFrame with merge_col as index
    """
    info_df = pd.read_csv(info_file, index_col=0, delimiter='\t')

    assert info_df.index.name == merge_col, \
        f'The file {info_file} is improperly formatted: The first column must be named "{merge_col}". ' \
        f'Current name: {info_df.index.name}. Remaining columns: {info_df.columns.tolist()}'

    if expected_overlap_set is not None:
        overlap_size = len(set.intersection(set(info_df.index), expected_overlap_set))
        if overlap_size == 0:
            logger.warning(f'The {merge_col}s in {info_file} do not match any {merge_col}s in {reference_file}')
        logger.debug(f'Loaded descriptions for {overlap_size} {merge_col}s')

    logger.debug(f'Loaded {merge_col} descriptions. columns={info_df.columns.tolist()}')
    assert not info_df.index.has_duplicates, \
        f'{info_file} contains duplicates: {info_df.index[info_df.index.duplicated()]}'
    return info_df


class MockCounter:
    """
    Imitate multiprocessing.Manager.Value / multiprocessing.managers.ValueProxy
    """

    def __init__(self):
        self._value: int = 0

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value


class MockLock:
    """
    Imitate multiprocessing.Manager.Lock / multiprocessing.managers.AcquirerProxy
    """

    def __init__(self):
        self._value = 0

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class AbstractNamespace:
    @classmethod
    def create_namespace(cls, ns, properties: {str: Any}):
        for name in cls.__dict__['__annotations__'].keys():
            setattr(ns, name, properties[name])
        return ns


def grasp_namespace(cls, ns):
    """
    This will copy the elements of the multiprocessing namespace into the "private" memory of the current process

    :param ns: multiprocessing.managers.Namespace
    :return: MockNameSpace
    """
    new_ns = cls()
    for name in cls.__dict__['__annotations__'].keys():
        value = getattr(ns, name)
        if name in ['lock', 'counter']:
            setattr(new_ns, name, value)
        else:
            setattr(new_ns, name, deepcopy(value))
    return new_ns


class AnalyzeTraitNamespace(AbstractNamespace):
    counter: MockCounter
    lock: MockLock
    outdir: str
    start_time: datetime
    genes_orig_df: pd.DataFrame
    genes_bool_df: pd.DataFrame
    gene_info_df: pd.DataFrame | None
    numeric_df: pd.DataFrame
    traits_df: pd.DataFrame
    trait_info_df: pd.DataFrame | None
    duplication_df: pd.DataFrame
    tree: object  #: ScoaryTree
    all_labels: set
    mt_f_method: str
    mt_f_cutoff: float
    max_genes: int
    worst_cutoff: None | float
    n_permut: int
    random_state: int
    pairwise: bool


class BinarizeTraitNamespace(AbstractNamespace):
    counter: MockCounter
    lock: MockLock
    outdir: str
    start_time: datetime
    numeric_df: pd.DataFrame
    random_state: int
    method: str
    alternative: str
    covariance_type: str
    cutoff: float
    random_state: int
