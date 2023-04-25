from functools import cache

from numba import njit
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from .ScoaryTree import ScoaryTree


def pick(
        tree: [],
        label_to_trait_a: {str: bool},
        trait_b_df: pd.DataFrame,
        calc_pvals: bool = True
) -> (np.array,):
    """
    Traverse the tree and perform pair picking

    :param tree: Tree in list form
    :param label_to_trait_a: maps each label of the tree to whether it has trait a
    :param trait_b_df: DataFrame (dtype:bool); columns: labels of the tree; rows: whether trait b is present
    :param calc_pvals: If False, binomial test will not be applied and best/worst will be None
    :return: (max_contr, max_suppo, max_oppos, best, worst) if calc_pvals else (max_contr, max_suppo, max_oppos)
    """

    assert not trait_b_df.isna().values.any()

    def _pick(left_label, right_label):
        # follow tree until terminal node
        if type(left_label) is not str:
            left = _pick(left_label[0], left_label[1])
        if type(right_label) is not str:
            right = _pick(right_label[0], right_label[1])

        # only load new leafs when needed for combination (safe RAM)
        if type(left_label) is str:
            left = init_leaf(
                trait_a=label_to_trait_a[left_label],
                trait_b_list=trait_b_df[left_label].to_numpy(dtype='bool')
            )
        if type(right_label) is str:
            right = init_leaf(
                trait_a=label_to_trait_a[right_label],
                trait_b_list=trait_b_df[right_label].to_numpy(dtype='bool')
            )

        combined = combine_branches(left, right)

        return combined

    values = _pick(tree[0], tree[1])

    max_contr = values[:, 0, :].max(axis=1)
    max_suppo = values[:, 1, :].max(axis=1)
    max_oppos = values[:, 2, :].max(axis=1)

    if not calc_pvals:
        return max_contr, max_suppo, max_oppos

    best, worst = apply_binomtest(max_contr, max_suppo, max_oppos)

    return max_contr, max_suppo, max_oppos, best, worst


def pick_single(
        tree: [],
        label_to_trait_a: {str: bool},
        label_to_trait_b: {str: bool},
        calc_pvals: bool = True
) -> {str: int | float}:
    res = pick(
        tree=tree,
        label_to_trait_a=label_to_trait_a,
        trait_b_df=pd.DataFrame([label_to_trait_b]),
        calc_pvals=calc_pvals
    )
    return dict(zip(
        ['max_contrasting_pairs', 'max_supporting_pairs', 'max_opposing_pairs', 'best_pval', 'worst_pval'],
        [v[0] for v in res]
    ))


def pick_nonrecursive(
        tree: [],
        label_to_trait_a: {str: bool},
        trait_b_df: pd.DataFrame,
        calc_pvals: bool = True
) -> (np.array, np.array, np.array, np.array, np.array):
    if tree.is_leaf:
        return init_leaf(
            trait_a=label_to_trait_a[tree.label],
            trait_b_list=trait_b_df[tree.label].to_numpy(dtype='bool')
        )

    stack = [[tree, 'right'], [tree, 'left']]

    while stack:
        current_parent, current_direction = stack[-1]
        current_node: ScoaryTree = getattr(current_parent, current_direction)

        if current_node.is_leaf:
            # current node is leaf
            this = init_leaf(
                trait_a=label_to_trait_a[current_node.label],
                trait_b_list=trait_b_df[current_node.label].to_numpy(dtype='bool')
            )

            # append data to parent
            current_node._values = this

            if current_direction == 'right':
                # found terminal node
                # # GO UP UNTIL CAN GO RIGHT
                while stack and stack[-1][1] == 'right':
                    ancestor_tree, ancestor_direction = stack.pop()
                    ancestor_tree._values = combine_branches(
                        ancestor_tree.left._values,
                        ancestor_tree.right._values
                    )
                    ancestor_tree.left._values = None
                    ancestor_tree.right._values = None

                if not stack:
                    # arrived at root node
                    break

            # pop left node -> go right next
            stack.pop()

        else:
            stack.extend([(current_node, 'right'), (current_node, 'left')])

    values = tree._values
    tree._values = None

    max_contr = values[:, 0, :].max(axis=1)
    max_suppo = values[:, 1, :].max(axis=1)
    max_oppos = values[:, 2, :].max(axis=1)

    if not calc_pvals:
        return max_contr, max_suppo, max_oppos

    best, worst = apply_binomtest(max_contr, max_suppo, max_oppos)

    return max_contr, max_suppo, max_oppos, best, worst


@cache
def _binomtest(k: int, n: int) -> float:
    # caching this function increases speed ~ 40x
    return binomtest(k=k, n=n).pvalue


def apply_binomtest(max_contr, max_suppo, max_oppos):
    n_traits = max_contr.shape[0]
    result = np.empty(shape=(2, n_traits), dtype='float')

    for i in range(n_traits):
        b = _binomtest(max_suppo[i], n=max_contr[i])
        w = _binomtest(max_oppos[i], n=max_contr[i])

        if b < w:
            result[0][i] = b
            result[1][i] = w
        else:
            result[0][i] = w
            result[1][i] = b
    return result


# selecting:values[<TRAITS>, <3 TYPES OF PAIRINGS>, <5 COMBINATIONS>]
# selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: 11; 1: 10; 2: 01; 3: 00; 4: nf>]

# values[n, 0, :] -> all max contrasting pairs for trait n
# values[n, 1, :] -> all max supporting pairs for trait n
# values[n, 2, :] -> all max opposing pairs for trait n

# values[n, 0, 0] -> max supporting pairs for trait n if condition '11' is added
# values[n, 0, 1] -> max supporting pairs for trait n if condition '10' is added
# values[n, 0, 2] -> max supporting pairs for trait n if condition '01' is added
# values[n, 0, 3] -> max supporting pairs for trait n if condition '00' is added
# values[n, 0, 4] -> max supporting pairs for trait n if condition 'nf' is added


@njit('int64[:, ::3, ::5](boolean, boolean[:])',
      cache=True, nogil=True, boundscheck=False, parallel=False)  # prange not better
def init_leaf(trait_a: bool, trait_b_list: np.array) -> np.array:
    n_traits = trait_b_list.shape[0]

    values = np.full(shape=(n_traits, 3, 5), fill_value=-1, dtype='int')
    if trait_a:
        values[:, :, 0][trait_b_list] = 0
        values[:, :, 1][~trait_b_list] = 0

    else:
        values[:, :, 2][trait_b_list] = 0
        values[:, :, 3][~trait_b_list] = 0

    return values


@njit('int64[::3, ::5], int64[::3, ::5]',
      cache=True, nogil=True, boundscheck=False, parallel=False)  # parallel kills performance
def calculate_max_nofree(left: np.array, right: np.array):
    values = np.full(shape=(3, 5), fill_value=-1, dtype='int')

    if left[0][4] > -1 and right[0][4] > -1:  # nf vs nf
        values[0][0] = left[0][4] + right[0][4]
        values[1][0] = left[1][4] + right[1][4]
        values[2][0] = left[2][4] + right[2][4]

    if left[0][0] > -1 and right[0][3] > -1:  # 11 vs 00
        values[0][1] = left[0][0] + right[0][3] + 1
        values[1][1] = left[1][0] + right[1][3] + 1
        values[2][1] = left[2][0] + right[2][3]

    if left[0][3] > -1 and right[0][0] > -1:  # 00 vs 11
        values[0][2] = left[0][3] + right[0][0] + 1
        values[1][2] = left[1][3] + right[1][0] + 1
        values[2][2] = left[2][3] + right[2][0]

    if left[0][1] > -1 and right[0][2] > -1:  # 10 vs 01
        values[0][3] = left[0][1] + right[0][2] + 1
        values[1][3] = left[1][1] + right[1][2]
        values[2][3] = left[2][1] + right[2][2] + 1

    if left[0][2] > -1 and right[0][1] > -1:  # 01 vs 10
        values[0][4] = left[0][2] + right[0][1] + 1
        values[1][4] = left[1][2] + right[1][1]
        values[2][4] = left[2][2] + right[2][1] + 1

    max_contr = values[0].max()

    max_suppo = -1
    for i in range(5):
        if values[0][i] == max_contr and values[1][i] > max_suppo:
            max_suppo = values[1][i]

    max_oppos = -1
    for i in range(5):
        if values[0][i] == max_contr and values[2][i] > max_oppos:
            max_oppos = values[2][i]

    return max_contr, max_suppo, max_oppos


@njit('int64, int64[::3, ::5], int64[::3, ::5]',
      cache=True, nogil=True, boundscheck=False, parallel=False)
def calculate_max_given_condition(condition: int, left: np.array, right: np.array):  # parallel kills performance
    values = np.full(shape=(3, 9), fill_value=-1, dtype='int')

    if left[0][condition] > -1:
        # compare condition with all conditions
        for i in range(5):
            values[0][i] = left[0][condition] + right[0][i]
            values[1][i] = left[1][condition] + right[1][i]
            values[2][i] = left[2][condition] + right[2][i]

    if right[0][condition] > -1:
        col_id = 5
        # compare all conditions with condition
        for i in range(5):
            if i == condition:  # this comparison has already been made above
                continue

            values[0][col_id] = left[0][i] + right[0][condition]
            values[1][col_id] = left[1][i] + right[1][condition]
            values[2][col_id] = left[2][i] + right[2][condition]

            col_id += 1

    max_contr = values[0].max()

    max_suppo = -1
    for i in range(9):
        if values[0][i] == max_contr and values[1][i] > max_suppo:
            max_suppo = values[1][i]

    max_oppos = -1
    for i in range(9):
        if values[0][i] == max_contr and values[2][i] > max_oppos:
            max_oppos = values[2][i]

    return max_contr, max_suppo, max_oppos


@njit('int64[:, ::3, ::5], int64[:, ::3, ::5]',
      cache=True, nogil=True, boundscheck=False, parallel=False)
def combine_branches(left: np.array, right: np.array):
    assert left.shape == right.shape
    n_traits = left.shape[0]

    values = np.full(shape=left.shape, fill_value=-1, dtype='int')

    # selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: 11; 1: 10; 2: 01; 3: 00; 4: nf>]
    for trait_id in range(n_traits):  # prange kills performance
        for cond in range(4):  # prange kills performance
            # {"11": 0, "10": 1, "01": 2, "00": 3, "nf": 4}
            max_contr, max_suppo, max_oppos = calculate_max_given_condition(
                cond,
                left[trait_id, :, :],
                right[trait_id, :, :]
            )
            values[trait_id, 0, cond] = max_contr
            values[trait_id, 1, cond] = max_suppo
            values[trait_id, 2, cond] = max_oppos
        max_contr, max_suppo, max_oppos = calculate_max_nofree(
            left[trait_id, :, :],
            right[trait_id, :, :]
        )
        values[trait_id, 0, 4] = max_contr
        values[trait_id, 1, 4] = max_suppo
        values[trait_id, 2, 4] = max_oppos

    return values
