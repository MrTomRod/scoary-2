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
        list_tree = scoary_tree.to_list

        self.assertEqual(expected_result, list_tree)

    def test_tree_from_genes_df(self):
        """
        Check if old scoary generates the equivalent tree based on genes presence/absence
        """
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
        # convert to ScoaryTree
        scoary_tree = ScoaryTree.from_presence_absence(genes_df)
        # convert to list
        list_tree = scoary_tree.to_list
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

    def test_uniquivy(self):
        label_to_trait = {'X': True, ' ': False}

        def apply(tree):
            return ScoaryTree.from_list(tree).uniquify(label_to_trait)

        expected_result = '(((01)1)(01))'
        for tree in (
                [['X', ' '], ['X', [' ', 'X']]],
                [[' ', 'X'], ['X', [' ', 'X']]],
                [[' ', 'X'], [[' ', 'X'], 'X']],
                [['X', [' ', 'X']], ['X', ' ']],
                [['X', ['X', ' ']], ['X', ' ']],
                [['X', ['X', ' ']], [' ', 'X']],
        ):
            unique_string = apply(tree)
            self.assertEqual(expected_result, unique_string)

        self.assertNotEqual(
            expected_result,
            apply([['X', ' '], ['X', [' ', ' ']]])
        )
