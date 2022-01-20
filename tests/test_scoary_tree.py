import time

from init_tests import *

from scoary.scoary import *

from scoary.ScoaryTree import *

from fast_fisher import odds_ratio

"""
todo:
- agreement between scoary 1
- max_supporting and best_picking do the same atm!
"""


class TestTreeFunctions(TestCase):
    def test_tree_from_list_to_list(self):
        expected_result = get_json('tetracycline', 'treelist')['as_list']
        # convert to ScoaryTree
        scoary_tree = ScoaryTree.from_list(expected_result)
        # convert back to list
        list_tree = scoary_tree.to_list()

        self.assertEqual(expected_result, list_tree)

    def test_tree_from_genes_df(self):
        """
        Check if old scoary generates the equivalent tree based on genes presence/absence
        """
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        # convert to ScoaryTree
        scoary_tree = ScoaryTree.from_presence_absence(genes_df)
        # convert to list
        list_tree = scoary_tree.to_list()
        # compare to Scoary 1
        expected_result = get_json('tetracycline', 'treelist')['as_list']
        self.assertTrue(is_equivalent_tree(expected_result, list_tree))

    def test_tree_from_newick_to_newick(self):
        """
        Check if newick tree is imported correctly
        """
        expected_result = get_json('tetracycline', 'treelist')['as_newick']
        scoary_tree = ScoaryTree.from_newick(newick=expected_result)
        newick = scoary_tree.to_newick()
        self.assertEqual(expected_result, newick)

    def test_pairs_simple(self):
        scoary_tree = ScoaryTree.from_list([['a', 'b'], ['c', 'd']])
        labels = scoary_tree.labels()
        assert labels == ['a', 'b', 'c', 'd']

        tree_data = {'a': (1, 1), 'b': (0, 0), 'c': (1, 1), 'd': (0, 0)}
        label_to_trait = {letter: trait for letter, (trait, gene) in tree_data.items()}
        label_to_gene = {letter: gene for letter, (trait, gene) in tree_data.items()}

        max_comparisons = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='contrasting')
        max_supporting = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='supporting')
        max_opposing = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='opposing')

        self.assertEqual(2, max_comparisons, msg='best picking of pairs failed')
        self.assertEqual(2, max_supporting, msg='best picking of pairs failed')
        self.assertEqual(0, max_opposing, msg='worst picking of pairs failed')

    def test_pairs_paper(self):
        scoary_tree = ScoaryTree.from_list(
            [[[[[[['1', '2'], ['3', '4']], '5'], '6'], '7'], '8'],
             [[[[['9', [['10', '11'], '12']], '13'], '14'], '15'], [['16', '17'], [['18', ['19', '20']], '21']]]]
        )
        print(scoary_tree)
        labels = scoary_tree.labels()
        assert labels == [str(v) for v in list(range(1, 22))]

        seq = [(0, 0), (0, 0), (1, 1), (1, 1), (1, 1), (0, 0), (0, 0), (1, 1), (1, 1), (0, 0), (1, 0), (0, 1), (0, 0), (1, 1), (1, 1), (0, 0), (0, 0),
               (1, 1), (1, 1), (0, 0), (1, 1), ]

        label_to_gene = {lab: bool(tup[0]) for tup, lab in zip(seq, labels)}
        label_to_trait = {lab: bool(tup[1]) for tup, lab in zip(seq, labels)}

        print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)

        max_comparisons = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='contrasting')
        max_supporting = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='supporting')
        max_opposing = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='opposing')

        # print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)

        self.assertEqual(7, max_supporting, msg='max_supporting of pairs failed')
        self.assertEqual(1, max_opposing, msg='max_opposing of pairs failed')
        self.assertEqual(7, max_comparisons, msg='max_comparisons of pairs failed')

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
            old_best_picking_pvalue = row.Best_pairwise_comp_p
            old_worst_picking_pvalue = row.Worst_pairwise_comp_p

            label_to_gene = genes_df.loc[gene].apply(bool).to_dict()

            max_comparisons = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='contrasting')
            max_supporting = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='supporting')
            max_opposing = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='opposing')

            print((old_max_comparisons, max_comparisons), (old_max_supporting, max_supporting), (old_max_opposing, max_opposing))

            if old_max_supporting != max_supporting or old_max_opposing != max_opposing:
                print(gene, scoary_tree)
                print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)
                assert False

    def test_scoary1_generated(self):
        ds_name = 'bigger_ds'  # 'small_ds'
        tree_id = 2
        genes_name = f'genes-{tree_id}'
        trait_name = f't2'
        trait_file = f't2-{tree_id}'

        genes_df = load_genes(get_path(ds_name, genes_name), delimiter=',', start_col=0)
        traits_df = load_traits(get_path(ds_name, 'traits'), delimiter=',')
        trait_df = traits_df[trait_name]

        expected_result = pd.read_csv(get_path(ds_name, trait_file))
        # new_odds_ratio = expected_result.apply(lambda row: odds_ratio(
        #     row['Number_pos_present_in'], row['Number_neg_present_in'], row['Number_pos_not_present_in'], row['Number_neg_not_present_in']), axis=1)
        # expected_result['Odds_ratio'] = new_odds_ratio
        # assert all(expected_result['Odds_ratio'] == new_odds_ratio)

        scoary_tree = ScoaryTree.from_presence_absence(genes_df)
        label_to_trait = traits_df.t1.apply(bool).to_dict()

        assert set(scoary_tree.labels()) == set(traits_df.index)

        for i, row in expected_result.sample(frac=1, random_state=42).iterrows():
            gene = row.Gene
            old_max_comparisons = row.Max_Pairwise_comparisons
            old_max_supporting = row.Max_supporting_pairs
            old_max_opposing = row.Max_opposing_pairs
            old_odds_ratio = row.Odds_ratio
            old_best_picking_pvalue = row.Best_pairwise_comp_p
            old_worst_picking_pvalue = row.Worst_pairwise_comp_p

            label_to_gene = genes_df.loc[gene].apply(bool).to_dict()
            if old_odds_ratio < 1:
                # only look at positive associations
                continue
                # invert gene depending on odds_ratio
                label_to_gene = {l: not g for l, g in label_to_gene.items()}
                inverted = True
            else:
                inverted = False

            max_comparisons = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='contrasting')
            max_supporting = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='supporting')
            max_opposing = count_max_pairings(scoary_tree, label_to_trait, label_to_gene, type='opposing')
            best_picking_pvalue = -1
            worst_picking_pvalue = -1

            print(gene, scoary_tree)
            print_tree_for_debugging(scoary_tree, label_to_gene, label_to_trait)

            df = pd.DataFrame({
                'Scoary': (old_max_comparisons, old_max_supporting, old_max_opposing,
                           # old_best_picking_pvalue, old_worst_picking_pvalue
                           ),
                'Scoary-2': (max_comparisons, max_supporting, max_opposing,
                             # best_picking_pvalue, worst_picking_pvalue
                             ),
            }, index=('max_comparisons', 'max_supporting', 'max_opposing',
                      # 'best_picking_pvalue', 'worst_picking_pvalue'
                      ), dtype=int)

            print(f'{inverted=} odds_ratio={old_odds_ratio} {gene=}')
            print(df)

            input()
            time.sleep(0.1)
