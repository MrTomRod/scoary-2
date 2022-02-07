from typing import Any, Callable
from unittest import TestCase

import os

# os.environ['PRINT_INIT'] = 'TRUE'
# os.environ['PRINT_COMBINE'] = 'TRUE'
# os.environ['PRINT_BEFORE_COMBINE'] = 'TRUE'
# os.environ['PRINT_AFTER_COMBINE'] = 'TRUE'

import numpy as np
import pandas as pd

from scoary.picking import pick

from scoary.scoary_1_picking import *

from timeit import default_timer as timer


def time_fn(fn: Callable, args=None, kwargs=None, n_times: int = 1) -> (float, Any):
    if kwargs is None:
        kwargs = {}
    if args is None:
        args = []

    diffs = []
    for i in range(n_times):
        start = timer()
        res = fn(*args, **kwargs)
        end = timer()
        diffs.append(end - start)  # Time in seconds, e.g. 5.38091952400282
    return np.mean(diffs), res


boolify = lambda t1, t2: f"{'A' if t1 else 'a'}{'B' if t2 else 'b'}"


def scoary_1_pick(tree: [], label_to_gene: {str: bool}, permuted_traits_df: pd.DataFrame):
    labels = set(permuted_traits_df.columns)

    max_contrasting = np.empty(shape=len(permuted_traits_df), dtype='int')
    max_supporting = np.empty(shape=len(permuted_traits_df), dtype='int')
    max_opposing = np.empty(shape=len(permuted_traits_df), dtype='int')

    for i, label_to_trait in permuted_traits_df.iterrows():
        gtc = {l: boolify(label_to_gene[l], label_to_trait[l]) for l in labels}
        phylo_tree, result_dict = convert_upgma_to_phylotree(tree, gtc)

        max_contrasting[i] = result_dict['Total']
        max_supporting[i] = result_dict['Pro']
        max_opposing[i] = result_dict['Anti']

    return max_contrasting, max_supporting, max_opposing


tree = [['strain1', 'strain2'], ['strain3', 'strain4']]

label_to_gene = {
    'strain1': True,
    'strain2': False,
    'strain3': False,
    'strain4': True,
}

permuted_traits_df = pd.DataFrame(
    [
        [True, True, False, False],
        [True, False, True, False],
        [True, False, False, True],
        [False, True, True, False],
        [False, True, False, True],
        [False, True, False, True],
        [False, True, False, True],
        [False, True, False, True],
    ], columns=['strain1', 'strain2', 'strain3', 'strain4']
)


class Test(TestCase):
    def test_old(self):
        for i in range(1000):
            mc, ms, mo = scoary_1_pick(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df)

    def test_new(self):
        for i in range(1000):
            mc, ms, mo = pick(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df)

    def test_identical(self):
        print('SCOARY 1')
        mc_1, ms_1, mo_1 = scoary_1_pick(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df)
        print('SCOARY 2')
        mc_2, ms_2, mo_2 = pick(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df)

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_pick(self):
        pick(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df)

        print('SCOARY 1')
        time_1, res = time_fn(
            scoary_1_pick,
            kwargs=dict(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df),
            n_times=1000
        )
        mc_1, ms_1, mo_1 = res

        print('SCOARY 2')
        time_2, res = time_fn(
            pick,
            kwargs=dict(tree=tree, label_to_gene=label_to_gene, permuted_traits_df=permuted_traits_df),
            n_times=1000
        )
        mc_2, ms_2, mo_2 = res

        print(time_1, time_2)
        print(time_1 / time_2, 'x improvement')

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')
