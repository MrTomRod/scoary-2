import time

import pandas as pd

from init_tests import *

from scoary.scoary import *

from scoary.ScoaryTree import *


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


    def test_prune(self):
        scoary_tree = ScoaryTree.from_list(
            [[[[[[['1', '2'], ['3', '4']], '5'], '6'], '7'], '8'],
             [[[[['9', [['10', '11'], '12']], '13'], '14'], '15'], [['16', '17'], [['18', ['19', '20']], '21']]]]
        )
        prune_labels = ['1', '2', '3', '18', '19', '21']
        pruned_tree = scoary_tree.prune(labels=prune_labels)
        real_labels = pruned_tree.labels()
        self.assertEqual(real_labels, prune_labels)

    def test_copy(self):
        scoary_tree = ScoaryTree.from_list(
            [[[[[[['1', '2'], ['3', '4']], '5'], '6'], '7'], '8'],
             [[[[['9', [['10', '11'], '12']], '13'], '14'], '15'], [['16', '17'], [['18', ['19', '20']], '21']]]]
        )
        copied_tree = scoary_tree.copy_nonrecursive()
        nonrec_copied_tree = scoary_tree.copy_nonrecursive()

        def confirm_copy(t1: ScoaryTree, t2: ScoaryTree):
            self.assertFalse(t1 is t2)
            if t1.is_leaf:
                self.assertTrue(t2.is_leaf)
                self.assertTrue(t1.label == t2.label)
            else:
                self.assertFalse(t2.is_leaf)
                confirm_copy(t1.left, t2.left)
                confirm_copy(t1.right, t2.right)

        confirm_copy(scoary_tree, copied_tree)
        confirm_copy(scoary_tree, nonrec_copied_tree)
        with self.assertRaises(AssertionError):
            confirm_copy(scoary_tree, scoary_tree)

    def test_prune_nonrecursive(self):
        scoary_tree = ScoaryTree.from_list(
            [[[[[[['1', '2'], ['3', '4']], '5'], '6'], '7'], '8'],
             [[[[['9', [['10', '11'], '12']], '13'], '14'], '15'], [['16', '17'], [['18', ['19', '20']], '21']]]]
        )
        prune_labels = ['1', '2', '3', '18', '19', '21']
        pruned_tree = scoary_tree.prune_nonrecursive(labels=prune_labels)
        real_labels = pruned_tree.labels()
        self.assertEqual(real_labels, prune_labels)


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
