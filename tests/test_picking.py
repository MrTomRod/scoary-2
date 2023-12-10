from typing import Any, Callable

from .init_tests import *

from scoary.utils import print_tree
from scoary.scoary import *
from scoary.analyze_trait import init_result_df, pair_picking
from scoary.ScoaryTree import ScoaryTree
from scoary.picking import pick, pick_nonrecursive, pick_single

from scoary.scoary_1_picking import *

from timeit import default_timer as timer

boolify = lambda t1, t2: f"{'A' if t1 else 'a'}{'B' if t2 else 'b'}"

dummy_tree = [['isolate1', 'isolate2'], ['isolate3', 'isolate4']]

dummy_trait_a = {
    'isolate1': True,
    'isolate2': False,
    'isolate3': False,
    'isolate4': True,
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
    ], columns=['isolate1', 'isolate2', 'isolate3', 'isolate4']
)


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


class Test(TestCase):
    def test_demo(self):
        tree = [['isolate1', 'isolate2'], [['isolate3', 'isolate4'], ['isolate5', 'isolate6']]]
        label_to_trait_a = {
            'isolate1': True,
            'isolate2': False,
            'isolate3': True,
            'isolate4': False,
            'isolate5': True,
            'isolate6': False,
        }
        label_to_trait_b = {
            'isolate1': True,
            'isolate2': False,
            'isolate3': True,
            'isolate4': False,
            'isolate5': True,
            'isolate6': False,
        }

        print_tree(
            ScoaryTree.from_list(tree),
            label_to_trait_a, label_to_trait_b
        )
        result = pick_single(tree, label_to_trait_a, label_to_trait_b, calc_pvals=True)
        print(result)

    def test_simple(self):
        mc_1, ms_1, mo_1 = scoary_1_pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df)
        mc_2, ms_2, mo_2 = pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df,
                                calc_pvals=False)

        self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
        self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
        self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_benchmark_tetracycline(self, run_scoary_1=True):
        # Scoary1 took 23.241052357999614 sec
        # Scoary2 took 0.49214521629996855 sec
        # Scoary1 vs Scoary2: 47.22397290118921x improvement
        tetr_tree = get_json('../data/tetracycline/expected_result.json')['as_list']
        _, tetr_genes_df = load_genes('../data/tetracycline/Gene_presence_absence.csv', gene_data_type='gene-count',
                                      ignore=roary_ignore)
        _, tetr_traits_df = load_traits('../data/tetracycline/Tetracycline_resistance.csv', trait_data_type='binary:,')

        tetr_label_to_gene = tetr_traits_df['Tetracycline_resistance'].to_dict()

        # jit compile
        pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df, calc_pvals=False)

        if run_scoary_1:
            print('Scoary1')
            time_1, res = time_fn(
                scoary_1_pick,
                kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df),
                n_times=5
            )
            mc_1, ms_1, mo_1 = res
        else:
            time_1 = 19.

        print('Scoary2')
        time_2, res = time_fn(
            pick,
            kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df,
                        calc_pvals=False),
            n_times=20
        )
        mc_2, ms_2, mo_2 = res

        print(f'Scoary1 took {time_1} sec')
        print(f'Scoary2 took {time_2} sec')
        print(f'Scoary1 vs Scoary2: {time_1 / time_2}x improvement')  # 33.88 x

        if run_scoary_1:
            self.assertTrue(all(np.equal(mc_1, mc_2)), msg='contrasting')
            self.assertTrue(all(np.equal(ms_1, ms_2)), msg='supporting')
            self.assertTrue(all(np.equal(mo_1, mo_2)), msg='opposing')

    def test_tetracycline_norecursive(self, run_scoary_1=True):
        # Scoary1 took 23.021255266200022 sec
        # Scoary2nonrec took 0.5782416850000118 sec
        # Scoary1 vs Scoary2nonrec: 39.81251415002976x improvement
        tetr_tree = get_json('../data/tetracycline/expected_result.json')['as_list']
        _, tetr_genes_df = load_genes('../data/tetracycline/Gene_presence_absence.csv', gene_data_type='gene-count',
                                      ignore=roary_ignore)
        _, tetr_traits_df = load_traits('../data/tetracycline/Tetracycline_resistance.csv', trait_data_type='binary:,')

        tetr_label_to_gene = tetr_traits_df['Tetracycline_resistance'].to_dict()

        # jit compile
        pick(tree=dummy_tree, label_to_trait_a=dummy_trait_a, trait_b_df=dummy_trait_b_df, calc_pvals=False)

        if run_scoary_1:
            print('Scoary1')
            time_1, res = time_fn(
                scoary_1_pick,
                kwargs=dict(tree=tetr_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df),
                n_times=5
            )
            mc_1, ms_1, mo_1 = res
        else:
            time_1 = 19.

        print('Scoary2')
        tetr_scoary_tree = ScoaryTree.from_list(tetr_tree)
        time_2, res = time_fn(
            pick_nonrecursive,
            kwargs=dict(tree=tetr_scoary_tree, label_to_trait_a=tetr_label_to_gene, trait_b_df=tetr_genes_df,
                        calc_pvals=False),
            n_times=20
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

        print_tree(scoary_tree, label_to_gene, label_to_trait)

        res = pick(
            scoary_tree.to_list,
            label_to_trait_a=label_to_trait,
            trait_b_df=pd.DataFrame(label_to_gene, index=['fakegene']),
            calc_pvals=False
        )

        max_comparisons = res[0][0]
        max_supporting = res[1][0]
        max_opposing = res[2][0]

        print_tree(scoary_tree, label_to_gene, label_to_trait)

        self.assertEqual(7, max_comparisons, msg='max_comparisons of pairs failed')
        self.assertEqual(7, max_supporting, msg='max_supporting of pairs failed')
        self.assertEqual(1, max_opposing, msg='max_opposing of pairs failed')

    def test_pairs_scoary1(self):
        _, genes_df = load_genes('../data/tetracycline/Gene_presence_absence.csv', gene_data_type='gene-count', ignore=roary_ignore)
        _, traits_df = load_traits('../data/tetracycline/Tetracycline_resistance.csv', trait_data_type='binary:,')
        expected_result = pd.read_csv('../data/tetracycline/fisher_permute100.results.csv')

        scoary_tree = ScoaryTree.from_presence_absence(genes_df)
        label_to_trait = traits_df.Tetracycline_resistance.apply(bool).to_dict()

        assert set(scoary_tree.labels()) == set(traits_df.index)
        assert not traits_df.Tetracycline_resistance.hasnans

        for i, row in expected_result.iterrows():
            gene = row.Gene
            print(gene)
            old_max_comparisons = row.Max_Pairwise_comparisons
            old_max_supporting = row.Max_supporting_pairs
            old_max_opposing = row.Max_opposing_pairs
            old_best = row.Best_pairwise_comp_p
            old_worst = row.Worst_pairwise_comp_p

            label_to_gene = genes_df.loc[gene].apply(bool).to_dict()

            res = pick(
                scoary_tree.to_list,
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
                    print_tree(scoary_tree, label_to_gene, label_to_trait)
                    self.fail(msg=f'Disagreement between Scoary1 and Scoary2')

    def test_scoary1_generated(self):
        _, genes_df = load_genes('../data/bigger_ds/pres_abs.csv', 'gene-count:,')
        _, traits_df = load_traits('../data/bigger_ds/trait_trees.csv', trait_data_type='binary:,')

        for trait_name in ['t1', 't2']:
            label_to_trait = traits_df[trait_name].apply(bool).to_dict()
            expected_result = pd.read_csv(f'../data/bigger_ds/{trait_name}.results.csv')

            with open('../data/bigger_ds/newick.nwk') as f:
                newick = f.read()
            scoary_tree = ScoaryTree.from_newick(newick)
            # scoary_tree = ScoaryTree.from_presence_absence(genes_df)

            result_df = init_result_df(genes_df, pd.Series(label_to_trait, dtype='boolean'))
            result_df = pair_picking(result_df, genes_df, scoary_tree, label_to_trait)

            assert set(scoary_tree.labels()) == set(traits_df.index)

            for i, row in expected_result.sample(frac=1, random_state=42).iterrows():
                gene = row.Gene
                old_max_comparisons = row.Max_Pairwise_comparisons
                old_max_supporting = row.Max_supporting_pairs
                old_max_opposing = row.Max_opposing_pairs
                old_best = row.Best_pairwise_comp_p
                old_worst = row.Worst_pairwise_comp_p

                new_row = result_df[result_df['Gene'] == gene].iloc[0]

                comparisons = {
                    'max_comparisons': (old_max_comparisons, new_row.contrasting),
                    'max_supporting': (old_max_supporting, new_row.supporting),
                    'max_opposing': (old_max_opposing, new_row.opposing),
                    'best': (old_best, new_row.best),
                    'worst': (old_worst, new_row.worst),
                }

                for comparison, (old, new) in comparisons.items():
                    if not np.isclose(old, new):
                        print(f'Error on {gene=} / {comparison=}')
                        print(comparisons)
                        label_to_gene = genes_df.loc[gene].apply(bool).to_dict()
                        print_tree(scoary_tree, label_to_gene, label_to_trait)
                        self.fail(msg=f'Disagreement between Scoary1 and Scoary2')
