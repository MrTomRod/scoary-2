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


def scoary_1_pick(tree: [], label_to_trait_a: {str: bool}, trait_b_df: pd.DataFrame):
    labels = set(trait_b_df.columns)

    max_contrasting = np.empty(shape=len(trait_b_df), dtype='int')
    max_supporting = np.empty(shape=len(trait_b_df), dtype='int')
    max_opposing = np.empty(shape=len(trait_b_df), dtype='int')

    for i, (_, label_to_trait) in enumerate(trait_b_df.iterrows()):
        gtc = {l: boolify(label_to_trait_a[l], label_to_trait[l]) for l in labels}
        phylo_tree, result_dict = convert_upgma_to_phylotree(tree, gtc)

        max_contrasting[i] = result_dict['Total']
        max_supporting[i] = result_dict['Pro']
        max_opposing[i] = result_dict['Anti']

    return max_contrasting, max_supporting, max_opposing


dummy_tree = [['strain1', 'strain2'], ['strain3', 'strain4']]

dummy_trait_a = {
    'strain1': True,
    'strain2': False,
    'strain3': False,
    'strain4': True,
}

dummy_trait_b_df = pd.DataFrame(
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
    def test_identical(self):
        print('SCOARY 1')
        mc_1, ms_1, mo_1 = scoary_1_pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df)
        print('SCOARY 2')
        mc_2, ms_2, mo_2 = pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df,
                                calc_pvals=False)

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_time(self):
        pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df, calc_pvals=False)

        print('SCOARY 1')
        time_1, res = time_fn(
            scoary_1_pick,
            kwargs=dict(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df),
            n_times=1000
        )
        mc_1, ms_1, mo_1 = res

        print('SCOARY 2')
        time_2, res = time_fn(
            pick,
            kwargs=dict(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df,
                        calc_pvals=False),
            n_times=1000
        )
        mc_2, ms_2, mo_2 = res

        print(time_1, time_2)
        print(time_1 / time_2, 'x improvement')

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_tetracycline(self, run_scoary_1=False):
        from scoary.scoary import load_genes, load_traits
        from tests.init_tests import get_json, get_path

        tetr_tree = get_json('tetracycline', 'treelist')['as_list']
        tetr_genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        tetr_traits_df = load_traits(get_path('tetracycline', 'traits'), delimiter=',')

        tetr_label_to_gene = tetr_traits_df['Tetracycline_resistance'].to_dict()

        # jit compile
        pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df, calc_pvals=False)

        if run_scoary_1:
            print('SCOARY 1')
            time_1, res = time_fn(
                scoary_1_pick,
                kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df),
                n_times=5
            )
            mc_1, ms_1, mo_1 = res
        else:
            time_1 = 19.

        print('SCOARY 2')
        time_2, res = time_fn(
            pick,
            kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df,
                        calc_pvals=False),
            n_times=20
        )
        mc_2, ms_2, mo_2 = res

        print(f'Scoary1 took {time_1} sec')
        print(f'Scoary2 took {time_2} sec')
        print(f'Scoary1 vs Scoary2: {time_1 / time_2}x improvement')

        if run_scoary_1:
            self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
            self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
            self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')
