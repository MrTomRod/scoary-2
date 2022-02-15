from .init_tests import *
from typing import Any, Callable

from scoary.scoary import load_genes, load_traits
from scoary.ScoaryTree import ScoaryTree
from scoary.picking import pick, pick_nonrecursive

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
        print('SCOARY 2 nonrecursive')
        mc_3, ms_3, mo_3 = pick_nonrecursive(tree=ScoaryTree.from_list(dummy_tree), label_to_trait_a=dummy_trait_a,
                                             trait_b_df=dummy_trait_b_df,
                                             calc_pvals=False)

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')
        self.assertTrue(all(np.equal(mc_1, mc_3)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_3)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_3)), msg='opposing')

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
            n_times=10
        )
        mc_2, ms_2, mo_2 = res

        print(f'Scoary1 took {time_1} sec')
        print(f'Scoary2 took {time_2} sec')
        print(f'Scoary1 vs Scoary2: {time_1 / time_2}x improvement')

        if run_scoary_1:
            self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
            self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
            self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_tetracycline_norecursive(self, run_scoary_1=False):
        from scoary.scoary import load_genes, load_traits
        from tests.init_tests import get_json, get_path

        tetr_tree = ScoaryTree.from_list(get_json('tetracycline', 'treelist')['as_list'])
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
            pick_nonrecursive,
            kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df,
                        calc_pvals=False),
            n_times=10
        )
        mc_2, ms_2, mo_2 = res

        print(f'Scoary1 took {time_1} sec')
        print(f'Scoary2nonrec took {time_2} sec')
        print(f'Scoary1 vs Scoary2nonrec: {time_1 / time_2}x improvement')

        if run_scoary_1:
            self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
            self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
            self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_pairs_paper(self):
        scoary_tree = ScoaryTree.from_list(
            [[[[[[['1', '2'], ['3', '4']], '5'], '6'], '7'], '8'],
             [[[[['9', [['10', '11'], '12']], '13'], '14'], '15'], [['16', '17'], [['18', ['19', '20']], '21']]]]
        )
        print(scoary_tree)
        labels = scoary_tree.labels()
        assert labels == [str(v) for v in list(range(1, 22))]

        seq = [(0, 0), (0, 0), (1, 1), (1, 1), (1, 1), (0, 0), (0, 0), (1, 1), (1, 1), (0, 0), (1, 0), (0, 1), (0, 0),
               (1, 1), (1, 1), (0, 0), (0, 0),
               (1, 1), (1, 1), (0, 0), (1, 1), ]

        label_to_gene = {lab: bool(tup[0]) for tup, lab in zip(seq, labels)}
        label_to_trait = {lab: bool(tup[1]) for tup, lab in zip(seq, labels)}

        print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)

        res = pick(
            scoary_tree.to_list(),
            label_to_trait_a=label_to_trait,
            trait_b_df=pd.DataFrame(label_to_gene, index=['fakegene']),
            calc_pvals=False
        )

        max_comparisons = res[0][0]
        max_supporting = res[1][0]
        max_opposing = res[2][0]

        print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)

        self.assertEqual(7, max_comparisons, msg='max_comparisons of pairs failed')
        self.assertEqual(7, max_supporting, msg='max_supporting of pairs failed')
        self.assertEqual(1, max_opposing, msg='max_opposing of pairs failed')

    def test_pairs_scoary1(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        traits_df = load_traits(get_path('tetracycline', 'traits'), delimiter=',')
        scoary1_res_df = pd.read_csv(get_path('tetracycline', 'scoary1-result'))

        scoary_tree = ScoaryTree.from_presence_absence(genes_df)
        label_to_trait = traits_df.Tetracycline_resistance.apply(bool).to_dict()

        assert set(scoary_tree.labels()) == set(traits_df.index)
        assert not traits_df.Tetracycline_resistance.hasnans

        for i, row in scoary1_res_df[10:].iterrows():
            gene = row.Gene
            old_max_comparisons = row.Max_Pairwise_comparisons
            old_max_supporting = row.Max_supporting_pairs
            old_max_opposing = row.Max_opposing_pairs
            old_odds_ratio = row.Odds_ratio
            old_best = row.Best_pairwise_comp_p
            old_worst = row.Worst_pairwise_comp_p

            label_to_gene = genes_df.loc[gene].apply(bool).to_dict()

            res = pick(
                scoary_tree.to_list(),
                label_to_trait_a=label_to_trait,
                trait_b_df=pd.DataFrame(label_to_gene, index=['fakegene']),
                calc_pvals=True
            )

            comparisons = {
                'max_comparisons': (old_max_comparisons, res[0][0]),
                'max_supporting': (old_max_supporting, res[1][0]),
                'max_opposing': (old_max_opposing, res[2][0]),
                'best': (old_best, res[3][0]),
                'worst': (old_worst, res[4][0]),
            }

            for comparison, (old, new) in comparisons.items():
                if not np.isclose(old, new):
                    print(gene, comparison, old, new, scoary_tree)
                    print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)
                    self.fail(msg=f'Disagreement between Scoary1 and Scoary2')
