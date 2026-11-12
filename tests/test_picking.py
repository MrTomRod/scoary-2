import numpy as np
import pandas as pd
import pytest
from scipy.stats import binomtest

from scoary.picking import pick, pick_single, pick_nonrecursive
from scoary.ScoaryTree import ScoaryTree
from scoary.scoary_1_picking import convert_upgma_to_phylotree


def scoary_1_pick(tree: list, label_to_trait_a: dict, trait_b_df: pd.DataFrame):
    """Reference implementation: the original Scoary 1 pair-picking algorithm."""
    boolify = lambda a, b: f"{'A' if a else 'a'}{'B' if b else 'b'}"
    res = []
    for _, label_to_trait_b in trait_b_df.iterrows():
        gtc = {l: boolify(label_to_trait_a[l], label_to_trait_b[l]) for l in trait_b_df.columns}
        _, d = convert_upgma_to_phylotree(tree, gtc)
        res.append((d['Total'], d['Pro'], d['Anti']))
    return tuple(np.array(x) for x in zip(*res))


def random_tree(labels: list, rng) -> list:
    labels = list(labels)
    rng.shuffle(labels)
    nodes = labels
    while len(nodes) > 1:
        i, j = sorted(rng.choice(len(nodes), 2, replace=False))
        b, a = nodes.pop(j), nodes.pop(i)
        nodes.append([a, b])
    return nodes[0]


TREE = [['i1', 'i2'], [['i3', 'i4'], ['i5', 'i6']]]
ALTERNATING = {'i1': True, 'i2': False, 'i3': True, 'i4': False, 'i5': True, 'i6': False}


def test_perfect_support():
    r = pick_single(TREE, ALTERNATING, ALTERNATING)
    assert (r['max_contrasting_pairs'], r['max_supporting_pairs'], r['max_opposing_pairs']) == (3, 3, 0)
    assert r['best_pval'] == pytest.approx(binomtest(3, 3).pvalue)
    assert r['worst_pval'] == pytest.approx(binomtest(3, 3).pvalue)


def test_perfect_opposition():
    inverse = {k: not v for k, v in ALTERNATING.items()}
    r = pick_single(TREE, ALTERNATING, inverse)
    assert (r['max_contrasting_pairs'], r['max_supporting_pairs'], r['max_opposing_pairs']) == (3, 0, 3)


def test_no_contrast_when_gene_is_constant():
    r = pick_single(TREE, ALTERNATING, {k: True for k in ALTERNATING}, calc_pvals=False)
    assert (r['max_contrasting_pairs'], r['max_supporting_pairs'], r['max_opposing_pairs']) == (0, 0, 0)


def test_calc_pvals_false_returns_three_arrays():
    res = pick(TREE, ALTERNATING, pd.DataFrame([ALTERNATING]), calc_pvals=False)
    assert len(res) == 3


def test_nan_is_rejected():
    df = pd.DataFrame([ALTERNATING], dtype=object)
    df.iloc[0, 0] = np.nan
    with pytest.raises(AssertionError):
        pick(TREE, ALTERNATING, df)


@pytest.mark.parametrize('seed', range(5))
def test_matches_scoary_1_and_nonrecursive(seed):
    rng = np.random.default_rng(seed)
    labels = [f'i{i}' for i in range(16)]
    tree = random_tree(labels, rng)
    trait_a = {l: bool(rng.integers(2)) for l in labels}
    trait_b_df = pd.DataFrame(rng.integers(0, 2, size=(30, len(labels))).astype(bool), columns=labels)

    expected = scoary_1_pick(tree, trait_a, trait_b_df)
    recursive = pick(tree, trait_a, trait_b_df, calc_pvals=False)
    nonrecursive = pick_nonrecursive(ScoaryTree.from_list(tree), trait_a, trait_b_df, calc_pvals=False)

    for name, e, r, n in zip(['contrasting', 'supporting', 'opposing'], expected, recursive, nonrecursive):
        np.testing.assert_array_equal(r, e, err_msg=f'pick vs scoary 1: {name}')
        np.testing.assert_array_equal(n, e, err_msg=f'pick_nonrecursive vs scoary 1: {name}')


def test_pvals_recursive_vs_nonrecursive():
    rng = np.random.default_rng(42)
    labels = [f'i{i}' for i in range(12)]
    tree = random_tree(labels, rng)
    trait_a = {l: bool(rng.integers(2)) for l in labels}
    trait_b_df = pd.DataFrame(rng.integers(0, 2, size=(10, len(labels))).astype(bool), columns=labels)
    r = pick(tree, trait_a, trait_b_df)
    n = pick_nonrecursive(ScoaryTree.from_list(tree), trait_a, trait_b_df)
    for a, b in zip(r, n):
        np.testing.assert_allclose(a, b)
