import pandas as pd

from init_tests import *

from unittest import TestCase
import numpy as np
from biotite.sequence.phylo import upgma as _biotite_upgma
from biotite.sequence.phylo.tree import TreeNode as BiotiteTreeNode

from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode

from scoary.upgma import upgma as scoary_upgma
from scoary.ScoaryTree import ScoaryTree


def biotite_upgma(tree: BiotiteTreeNode, labels: [str]) -> ScoaryTree:
    def convert(node: BiotiteTreeNode) -> ScoaryTree:
        """recursive function"""
        if node.is_leaf():
            return ScoaryTree(label=str(node_to_label[node]))
        else:
            return ScoaryTree(left=convert(node.children[0]), right=convert(node.children[1]))

    node_to_label: {BiotiteTreeNode: str} = {node: label for node, label in zip(tree.leaves, labels)}
    return convert(tree.root)


def scipy_upgma(distances, labels: [str]):
    """
    scipy.cluster.hierarchy.linkage: method=’average’ is called the UPGMA algorithm!
    """

    def convert(node: ClusterNode) -> ScoaryTree:
        """recursive function"""
        if node.count == 1:  # is_leaf
            return ScoaryTree(label=str(node_to_label[node.id]))
        else:
            return ScoaryTree(left=convert(node.left), right=convert(node.right))

    node_to_label: {BiotiteTreeNode: str} = dict(enumerate(labels))

    Z = linkage(distances, method='average')
    tree = to_tree(Z, False)
    return convert(tree)


class Test(TestCase):
    def test_upgma(self):
        distances = np.array([
            [0, 1, 7, 7, 9],
            [1, 0, 7, 6, 8],
            [7, 7, 0, 2, 4],
            [7, 6, 2, 0, 3],
            [9, 8, 4, 3, 0],
        ])
        labels = [f'l{i}' for i in range(5)]

        _biotite_tree = _biotite_upgma(distances)
        biotite_tree = biotite_upgma(_biotite_tree, labels=labels).to_list

        scipy_tree = scipy_upgma(distances, labels).to_list

        distances_df = pd.DataFrame(distances, columns=labels)
        scoary_tree = scoary_upgma(distances_df)

        print(biotite_tree)
        print(scoary_tree)
        print(scipy_tree)

        assert is_equivalent_tree(biotite_tree, scipy_tree)
        assert is_equivalent_tree(biotite_tree, scoary_tree)

    def test_scoary(self, size=20, n_tests=1000):
        labels = [f'l{i}' for i in range(size)]

        n_failures = 0
        for i in range(n_tests):
            matrix = np.random.randint(0, 2000, size=(size, size))
            symmetrical_matrix = (matrix + matrix.T) / 2

            distances_df = pd.DataFrame(symmetrical_matrix, columns=labels)
            scoary_tree = scoary_upgma(distances_df)

            _biotite_tree = _biotite_upgma(symmetrical_matrix)
            biotite_tree = biotite_upgma(_biotite_tree, labels=labels).to_list

            if not is_equivalent_tree(biotite_tree, scoary_tree):
                print('no match:')
                print(f'  {biotite_tree=}')
                print(f'   {scoary_tree=}')
                n_failures += 1

        print(f'{n_failures=} out of {n_tests} tests')
        print(n_failures / n_tests)
        self.assertLess(n_failures / n_tests, 0.05, f'Lots of failures, wtf?')

    def test_scipy(self, size=20, n_tests=1000):
        """
        0 % agreement!
        Evidently, scipy's upgma implementation works rather differently!
        """
        labels = [f'l{i}' for i in range(size)]

        n_failures = 0
        for i in range(n_tests):
            matrix = np.random.randint(0, 2000, size=(size, size))
            symmetrical_matrix = (matrix + matrix.T) / 2

            _biotite_tree = _biotite_upgma(symmetrical_matrix)
            biotite_tree = biotite_upgma(_biotite_tree, labels=labels).to_list

            scipy_tree = scipy_upgma(symmetrical_matrix, labels).to_list

            if not is_equivalent_tree(biotite_tree, scipy_tree):
                print('no match:')
                print(f'  {biotite_tree=}')
                print(f'    {scipy_tree=}')
                n_failures += 1

        print(f'{n_failures=} out of {n_tests} tests')
        print(n_failures / n_tests)
        self.assertLess(n_failures / n_tests, 0.05, f'Lots of failures, wtf?')
