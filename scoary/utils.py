import os
import sys
import logging
import warnings
from typing import Type
from datetime import datetime
import pandas as pd

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

ALLOWED_CORRECTIONS = {'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by',
                       'fdr_tsbh', 'fdr_tsbky'}


def setup_outdir(outdir: str) -> str:
    outdir = outdir.rstrip('/')
    assert not os.path.exists(outdir), f'ERROR: {outdir=} already exists!'
    os.makedirs(f'{outdir}/traits')
    return outdir


def setup_logging(path: str):
    # print only warnings to console
    # todo: is this needed?
    # stdout = logging.StreamHandler()
    # stdout.setLevel(logging.WARNING)
    # logging.getLogger('root').addHandler(stdout)

    # create logfile that contains all logs
    logfile = logging.FileHandler(path)
    logfile.setLevel(logging.INFO)

    # add handler to all existing loggers
    for name in logging.root.manager.loggerDict:
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.addHandler(logfile)


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
        logging.info(f'Setting new recursion limit: {self.old} -> {self.new}')
        sys.setrecursionlimit(self.new)

    def __exit__(self, *args, **kwargs):
        logging.info(f'Setting old recursion limit: {self.new} -> {self.old}')
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


def get_label_to_trait(trait: pd.Series) -> {str: bool}:
    return {l: bool(t) for l, t in trait.items() if not pd.isna(t)}


def get_all_label_to_gene(genes_df: pd.DataFrame) -> pd.DataFrame:
    assert not genes_df.isna().values.any(), 'genes_df contains NaN'
    return genes_df.T.applymap(bool).to_dict()


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
        f'The file {info_file} is improperly formatted: The first column must be named "{merge_col}".' \
        f'Current name: {info_df.index.name}. Remaining columns: {info_df.columns.tolist()}'

    if expected_overlap_set is not None:
        overlap_size = len(set.intersection(set(info_df.index), expected_overlap_set))
        assert overlap_size > 0, f'The {merge_col}s in {info_file} do not match any {merge_col}s in {reference_file}'
        logger.info(f'Loaded descriptions for {overlap_size} {merge_col}s')

    logger.info(f'Loaded {merge_col} descriptions. columns={info_df.columns.tolist()}')
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


class MockNamespace:
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
    mt_p_method: str
    mt_p_cutoff: float
    n_permut: int
    all_label_to_gene: dict
    random_state: int
    no_pairwise: bool


def create_namespace(ns, counter, lock, outdir,
                     genes_orig_df, genes_bool_df, gene_info_df,
                     numeric_df, traits_df, trait_info_df,
                     duplication_df, tree, all_labels,
                     mt_f_method, mt_f_cutoff, mt_p_method, mt_p_cutoff,
                     n_permut, all_label_to_gene, random_state, no_pairwise):
    ns.counter = counter
    ns.lock = lock
    ns.outdir = outdir
    ns.start_time = datetime.now()
    ns.genes_orig_df = genes_orig_df
    ns.genes_bool_df = genes_bool_df
    ns.gene_info_df = gene_info_df
    ns.numeric_df = numeric_df
    ns.traits_df = traits_df
    ns.trait_info_df = trait_info_df
    ns.duplication_df = duplication_df
    ns.tree = tree
    ns.all_labels = all_labels
    ns.mt_f_method = mt_f_method
    ns.mt_f_cutoff = mt_f_cutoff
    ns.mt_p_method = mt_p_method
    ns.mt_p_cutoff = mt_p_cutoff
    ns.n_permut = n_permut
    ns.all_label_to_gene = all_label_to_gene
    ns.random_state = random_state
    ns.no_pairwise = no_pairwise
    return ns
