import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

from scoary.upgma import upgma
from conftest import is_equivalent_tree


def scipy_upgma(distances: np.ndarray, labels: list) -> list:
    """Reference: scipy's 'average' linkage is UPGMA."""

    def convert(node):
        return labels[node.id] if node.is_leaf() else [convert(node.left), convert(node.right)]

    return convert(to_tree(linkage(squareform(distances, checks=False), method='average')))


def test_textbook_example():
    distances = np.array([
        [0, 1, 7, 7, 9],
        [1, 0, 7, 6, 8],
        [7, 7, 0, 2, 4],
        [7, 6, 2, 0, 3],
        [9, 8, 4, 3, 0],
    ])
    labels = [f'l{i}' for i in range(5)]
    tree = upgma(pd.DataFrame(distances, columns=labels))
    assert is_equivalent_tree(tree, [['l0', 'l1'], [['l2', 'l3'], 'l4']])
    assert is_equivalent_tree(tree, scipy_upgma(distances, labels))


@pytest.mark.parametrize('seed', range(20))
def test_matches_scipy_on_random_matrices(seed):
    rng = np.random.default_rng(seed)
    size = 15
    # continuous random distances -> no ties, so the topology is unique
    points = rng.random((size, 4))
    distances = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(-1))
    labels = [f'l{i}' for i in range(size)]
    tree = upgma(pd.DataFrame(distances, columns=labels))
    assert is_equivalent_tree(tree, scipy_upgma(distances, labels))
