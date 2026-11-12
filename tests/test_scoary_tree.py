import pandas as pd
import pytest

from scoary.ScoaryTree import ScoaryTree
from conftest import is_equivalent_tree

TREE = [['A', 'B'], ['C', ['D', 'E']]]


def test_list_round_trip():
    assert ScoaryTree.from_list(TREE).to_list == TREE


def test_newick_round_trip():
    tree = ScoaryTree.from_list(TREE)
    newick = tree.to_newick()
    assert newick == '((A,B),(C,(D,E)));'
    assert ScoaryTree.from_newick(newick).to_list == TREE


def test_labels():
    assert ScoaryTree.from_list(TREE).labels() == ['A', 'B', 'C', 'D', 'E']


def test_leaf():
    leaf = ScoaryTree(label='A')
    assert leaf.is_leaf
    assert leaf.labels() == ['A']


def test_copy_is_deep():
    tree = ScoaryTree.from_list(TREE)
    copy = tree.copy()
    assert copy.to_list == tree.to_list
    assert copy is not tree
    assert copy.left is not tree.left
    assert copy.right.right.left is not tree.right.right.left


@pytest.mark.parametrize('keep, expected', [
    (['A', 'B', 'C', 'D', 'E'], TREE),
    (['A', 'C', 'D', 'E'], ['A', ['C', ['D', 'E']]]),
    (['A', 'B', 'D'], [['A', 'B'], 'D']),
    (['B', 'E'], ['B', 'E']),
])
def test_prune(keep, expected):
    tree = ScoaryTree.from_list(TREE)
    pruned = tree.prune(keep)
    assert pruned.to_list == expected
    assert tree.to_list == TREE, 'pruning must not modify the original tree'


def test_prune_unknown_label_fails():
    with pytest.raises(AssertionError):
        ScoaryTree.from_list(TREE).prune(['A', 'B', 'not-in-tree'])


def test_rename():
    tree = ScoaryTree.from_list(TREE)
    assert tree.rename(str.lower).labels() == ['a', 'b', 'c', 'd', 'e']
    assert tree.labels() == ['A', 'B', 'C', 'D', 'E'], 'renaming must not modify the original tree'


def test_uniquify_is_order_independent():
    traits = {'A': True, 'B': False, 'C': True, 'D': False, 'E': True}
    t1 = ScoaryTree.from_list([['A', 'B'], ['C', ['D', 'E']]])
    t2 = ScoaryTree.from_list([[['E', 'D'], 'C'], ['B', 'A']])
    assert t1.uniquify(traits) == t2.uniquify(traits)


def test_from_presence_absence_groups_identical_isolates():
    # A/B and C/D have identical gene content, E is an outlier
    genes = pd.DataFrame({
        'A': [1, 1, 0, 0, 1],
        'B': [1, 1, 0, 0, 1],
        'C': [0, 0, 1, 1, 1],
        'D': [0, 0, 1, 1, 1],
        'E': [1, 0, 1, 0, 0],
    }).astype(bool)
    tree = ScoaryTree.from_presence_absence(genes)
    assert sorted(tree.labels()) == list('ABCDE')

    def find_cherry(t, a, b):
        if t.is_leaf:
            return False
        if t.left.is_leaf and t.right.is_leaf and {t.left.label, t.right.label} == {a, b}:
            return True
        return find_cherry(t.left, a, b) or find_cherry(t.right, a, b)

    assert find_cherry(tree, 'A', 'B')
    assert find_cherry(tree, 'C', 'D')


def test_tetracycline_tree(tetracycline_tree_list, tetracycline_genes):
    tree = ScoaryTree.from_list(tetracycline_tree_list)
    assert set(tree.labels()) == set(tetracycline_genes.columns)
    assert is_equivalent_tree(ScoaryTree.from_newick(tree.to_newick()).to_list, tetracycline_tree_list)
