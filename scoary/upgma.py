import numpy as np
import pandas as pd


def _find_min(arr: np.ndarray, n_cols: int) -> (int, int):
    min_index = int(np.nanargmin(arr))
    x, y = (min_index // n_cols, min_index % n_cols)
    assert x > y
    return x, y


def _merge(arr: np.ndarray, node_list: [], cluster_sizes: [int], x: int, y: int):
    n_rows, n_cols = arr.shape
    assert n_rows == n_cols == len(node_list) == len(cluster_sizes)
    assert x > y

    # update arr
    for p in range(len(arr)):
        if p in (x, y):
            continue
        px1, px2 = (x, p) if p < x else (p, x)
        py1, py2 = (y, p) if p < y else (p, y)
        assert not (np.isnan(arr[px1, px2]) or np.isnan(arr[py1, py2]))

        # calculate mean difference
        arr[py1, py2] = (arr[px1, px2] * cluster_sizes[x] + arr[py1, py2] * cluster_sizes[y]) / (cluster_sizes[x] + cluster_sizes[y])  # row, col

    # remove row and col x
    arr = np.delete(arr, x, 0)
    arr = np.delete(arr, x, 1)

    # update labels
    new_label = [node_list[y], node_list[x]]
    del node_list[x]
    node_list[y] = new_label

    # update cluster_sizes
    cluster_sizes[y] = cluster_sizes[x] + cluster_sizes[y]
    del cluster_sizes[x]

    assert arr.shape == (n_rows - 1, n_cols - 1)
    assert len(node_list) == n_rows - 1
    return arr, node_list, cluster_sizes


def _upgma(arr: np.ndarray, node_list: [str]) -> []:
    # fill triu with nan
    arr[np.triu_indices(arr.shape[0], k=0)] = np.nan

    cluster_sizes = [1 for _ in range(len(node_list))]

    while len(node_list) > 1:
        n_rows, n_cols = arr.shape
        assert len(node_list) == n_rows == n_cols

        # find next columns to merge
        x, y = _find_min(arr, len(node_list))

        # merge columns and update labels
        arr, node_list, cluster_sizes = _merge(arr, node_list, cluster_sizes, x, y)

    assert len(node_list) == 1
    tree = node_list[0]

    return tree


def upgma(distances: pd.DataFrame) -> []:
    """
    Apply UPGMA (unweighted pair group method with arithmetic mean) algorithm.
    Returns unweighted tree in nested list form.

    Insipred by 'Creating a Phylogenetic Tree' by Oxford Academic (https://www.youtube.com/watch?v=09eD4A_HxVQ)

    :param distances: pandas.DataFrame: values: symmetric array; columns: tree labels
    :return: tree: nested list of strings
    """
    # split pandas.DataFrame into numpy.ndarray and labels
    labels: [] = [str(c) for c in distances.columns]
    arr: np.ndarray = distances.values.astype(float)

    # sanity checks
    assert len(set(labels)) == len(labels), f'labels are not unique! {labels=}'
    assert np.allclose(arr, arr.T, rtol=1e-05, atol=1e-08), f'arr is not symmetric! arr:\n{distances.to_string()}'
    assert arr.shape[0] == arr.shape[1]
    assert not np.isnan(arr).any(), 'Distance matrix contains nan'
    assert not np.isinf(arr).any(), 'Distance matrix contains inf'
    assert not np.any(distances < 0), 'Distances must be positive'

    return _upgma(arr=arr, node_list=labels)
