from functools import cache

from numba import njit, prange
import numpy as np
import pandas as pd
from scipy.stats import binom_test

from .ScoaryTree import ScoaryTree


def pick(
        tree: [],
        label_to_trait_a: {str: bool},
        trait_b_df: pd.DataFrame,
        calc_pvals: bool = True
) -> (np.array, np.array, np.array, np.array, np.array):
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
        if type(left_label) is str:
            assert type(label_to_trait_a[left_label]) is bool
            left = init_leaf(
                trait_a=label_to_trait_a[left_label],
                trait_b_list=trait_b_df[left_label].to_numpy(dtype='bool')
            )
        else:
            left = _pick(left_label[0], left_label[1])
        if type(right_label) is str:
            right = init_leaf(
                trait_a=label_to_trait_a[right_label],
                trait_b_list=trait_b_df[right_label].to_numpy(dtype='bool')
            )
        else:
            right = _pick(right_label[0], right_label[1])

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
def _binomtest(x: int, n: int) -> float:
    # caching this function increases speed ~ 40x
    # older (deprecated) function binom_test is faster than newer binomtest
    return binom_test(x, n=n)


def apply_binomtest(max_contr, max_suppo, max_oppos):
    n_traits = max_contr.shape[0]
    result = np.empty(shape=(2, n_traits), dtype='float')

    for i in prange(n_traits):
        b = _binomtest(max_suppo[i], n=max_contr[i])
        w = _binomtest(max_contr[i] - max_oppos[i], n=max_contr[i])

        if b < w:
            result[0][i] = b
            result[1][i] = w
        else:
            result[0][i] = w
            result[1][i] = b
    return result


# selecting:values[<TRAITS>, <3 TYPES OF PAIRINGS>, <5 COMBINATIONS>]
# selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: 11; 1: 10; 2: 01; 3: 00; 4: free>]

# values[n, 0, :] -> all max contrasting pairs for trait n
# values[n, 1, :] -> all max supporting pairs for trait n
# values[n, 2, :] -> all max opposing pairs for trait n

# values[n, 0, 0] -> max supporting pairs for trait n if condition '11' is added
# values[n, 0, 1] -> max supporting pairs for trait n if condition '10' is added
# values[n, 0, 2] -> max supporting pairs for trait n if condition '01' is added
# values[n, 0, 3] -> max supporting pairs for trait n if condition '00' is added
# values[n, 0, 4] -> max supporting pairs for trait n if condition 'free' is added


@njit('int64[:, ::3, ::5](b1, boolean[:])', nogil=True, boundscheck=False, parallel=False)  # prange not better
def init_leaf(trait_a: bool, trait_b_list: np.array) -> np.array:
    n_traits = trait_b_list.shape[0]

    values = np.full(shape=(n_traits, 3, 5), fill_value=-1, dtype='int')
    if trait_a:
        for i, trait in enumerate(trait_b_list):
            if trait:
                values[i, 0, 0] = 0  # contrasting pairs
                values[i, 1, 0] = 0  # supporting pairs
                values[i, 2, 0] = 0  # opposing pairs
            else:
                values[i, 0, 1] = 0  # pairs_contr
                values[i, 1, 1] = 0  # pairs_supporting
                values[i, 2, 1] = 0  # opposing pairs
    else:
        for i, trait in enumerate(trait_b_list):
            if trait:
                values[i, 0, 2] = 0  # contrasting pairs
                values[i, 1, 2] = 0  # supporting pairs
                values[i, 2, 2] = 0  # opposing pairs
            else:
                values[i, 0, 3] = 0  # contrasting pairs
                values[i, 1, 3] = 0  # supporting pairs
                values[i, 2, 3] = 0  # opposing pairs

    return values


@njit('int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False, parallel=False)  # parallel kills performance
def calculate_max_free(left: np.array, right: np.array):
    values = np.full(shape=(3, 5), fill_value=-1, dtype='int')

    if left[0][4] > -1 and right[0][4] > -1:  # 0 vs 0
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
    for i in prange(5):
        if values[0][i] == max_contr and values[1][i] > max_suppo:
            max_suppo = values[1][i]

    max_oppos = -1
    for i in prange(5):
        if values[0][i] == max_contr and values[2][i] > max_oppos:
            max_oppos = values[2][i]

    return max_contr, max_suppo, max_oppos


@njit('int64, int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False, parallel=False)
def calculate_max_given_condition(condition: int, left: np.array, right: np.array):  # parallel kills performance
    values = np.full(shape=(3, 9), fill_value=-1, dtype='int')

    if left[0][condition] > -1:
        # compare condition with all conditions
        for i in prange(5):
            values[0][i] = left[0][condition] + right[0][i]
            values[1][i] = left[1][condition] + right[1][i]
            values[2][i] = left[2][condition] + right[2][i]

    if right[0][condition] > -1:
        col_id = 5
        # compare all conditions with condition
        for i in prange(5):
            if i == condition:  # this comparison has already been made 00ove
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


@njit('int64[:, ::3, ::5], int64[:, ::3, ::5]', nogil=True, boundscheck=False, parallel=False)
def combine_branches(left: np.array, right: np.array):
    assert left.shape == right.shape
    n_traits = left.shape[0]

    values = np.full(shape=left.shape, fill_value=-1, dtype='int')

    # selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: 11; 1: 10; 2: 01; 3: 00; 4: free>]
    for trait_id in prange(n_traits):  # prange kills performance
        for cond in prange(4):  # prange kills performance
            # {"11": 0, "10": 1, "01": 2, "00": 3, "0": 4}
            max_contr, max_suppo, max_oppos = calculate_max_given_condition(
                cond,
                left[trait_id, :, :],
                right[trait_id, :, :]
            )
            values[trait_id, 0, cond] = max_contr
            values[trait_id, 1, cond] = max_suppo
            values[trait_id, 2, cond] = max_oppos
        max_contr, max_suppo, max_oppos = calculate_max_free(
            left[trait_id, :, :],
            right[trait_id, :, :]
        )
        values[trait_id, 0, 4] = max_contr
        values[trait_id, 1, 4] = max_suppo
        values[trait_id, 2, 4] = max_oppos

    return values

# @njit('int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False)  # no prange
# def calculate_max_free_SC1(left: np.array, right: np.array):
#     max_contr_free = -1
#     max_contr_1100 = -1
#     max_contr_0011 = -1
#     max_contr_1001 = -1
#     max_contr_0110 = -1
#
#     max_suppo_free = -1
#     max_suppo_1100 = -1
#     max_suppo_0011 = -1
#     max_suppo_1001 = -1
#     max_suppo_0110 = -1
#
#     max_oppos_free = -1
#     max_oppos_1100 = -1
#     max_oppos_0011 = -1
#     max_oppos_1001 = -1
#     max_oppos_0110 = -1
#
#     l_max_contr = left[0]
#     r_max_contr = right[0]
#
#     l_max_suppo = left[1]
#     r_max_suppo = right[1]
#
#     l_max_oppos = left[2]
#     r_max_oppos = right[2]
#
#     if l_max_contr[4] > -1 and r_max_contr[4] > -1:
#         max_contr_free = l_max_contr[4] + r_max_contr[4]
#         max_suppo_free = l_max_suppo[4] + r_max_suppo[4]
#         max_oppos_free = l_max_oppos[4] + r_max_oppos[4]
#
#     if l_max_contr[0] > -1 and r_max_contr[3] > -1:
#         max_contr_1100 = l_max_contr[0] + r_max_contr[3] + 1
#         max_suppo_1100 = l_max_suppo[0] + r_max_suppo[3] + 1
#         max_oppos_1100 = l_max_oppos[0] + r_max_oppos[3]
#
#     if l_max_contr[3] > -1 and r_max_contr[0] > -1:
#         max_contr_0011 = l_max_contr[3] + r_max_contr[0] + 1
#         max_suppo_0011 = l_max_suppo[3] + r_max_suppo[0] + 1
#         max_oppos_0011 = l_max_oppos[3] + r_max_oppos[0]
#
#     if l_max_contr[1] > -1 and r_max_contr[2] > -1:
#         max_contr_1001 = l_max_contr[1] + r_max_contr[2] + 1
#         max_suppo_1001 = l_max_suppo[1] + r_max_suppo[2]
#         max_oppos_1001 = l_max_oppos[1] + r_max_oppos[2] + 1
#
#     if l_max_contr[2] > -1 and r_max_contr[1] > -1:
#         max_contr_0110 = l_max_contr[2] + r_max_contr[1] + 1
#         max_suppo_0110 = l_max_suppo[2] + r_max_suppo[1]
#         max_oppos_0110 = l_max_oppos[2] + r_max_oppos[1] + 1
#
#     max_contr = max(max_contr_free, max_contr_1100,
#                     max_contr_0011, max_contr_1001, max_contr_0110)
#
#     # Calculate max number of propairs
#     max_suppo = -1  # max_suppo can never go below -1
#     if max_contr == max_contr_free:
#         max_suppo = max(max_suppo, max_suppo_free)
#     if max_contr == max_contr_1100:
#         max_suppo = max(max_suppo, max_suppo_1100)
#     if max_contr == max_contr_0011:
#         max_suppo = max(max_suppo, max_suppo_0011)
#     if max_contr == max_contr_1001:
#         max_suppo = max(max_suppo, max_suppo_1001)
#     if max_contr == max_contr_0110:
#         max_suppo = max(max_suppo, max_suppo_0110)
#
#     # Calculate max number of antipairs
#     max_oppos = -1  # max_oppos can never go below -1
#     if max_contr == max_contr_free:
#         max_oppos = max(max_oppos, max_oppos_free)
#     if max_contr == max_contr_1100:
#         max_oppos = max(max_oppos, max_oppos_1100)
#     if max_contr == max_contr_0011:
#         max_oppos = max(max_oppos, max_oppos_0011)
#     if max_contr == max_contr_1001:
#         max_oppos = max(max_oppos, max_oppos_1001)
#     if max_contr == max_contr_0110:
#         max_oppos = max(max_oppos, max_oppos_0110)
#
#     return max_contr, max_suppo, max_oppos
#
#
# @njit('int64, int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False, parallel=False)
# def calculate_max_given_condition_SC1(condition, left: np.array, right: np.array):  # calculate_max_condition
#     other_conditions = [c for c in range(4) if c != condition]
#
#     max_contr_1, max_contr_2, max_contr_3, max_contr_4, \
#     max_contr_5, max_contr_6, max_contr_7, max_contr_8, max_contr_9 = [-1] * 9
#
#     max_suppo_1, max_suppo_2, max_suppo_3, max_suppo_4, \
#     max_suppo_5, max_suppo_6, max_suppo_7, max_suppo_8, max_suppo_9 = [-1] * 9
#
#     max_oppos_1, max_oppos_2, max_oppos_3, max_oppos_4, \
#     max_oppos_5, max_oppos_6, max_oppos_7, max_oppos_8, max_oppos_9 = [-1] * 9
#
#     l_max_contr = left[0]
#     r_max_contr = right[0]
#
#     l_max_suppo = left[1]
#     r_max_suppo = right[1]
#
#     l_max_oppos = left[2]
#     r_max_oppos = right[2]
#
#     # special case: comparison with nonfree
#     if l_max_contr[condition] > -1 and r_max_contr[4] > -1:
#         max_contr_1 = l_max_contr[condition] + r_max_contr[4]
#         max_suppo_1 = l_max_suppo[condition] + r_max_suppo[4]
#         max_oppos_1 = l_max_oppos[condition] + r_max_oppos[4]
#
#     # comparison with all others
#     if l_max_contr[condition] > -1 and r_max_contr[other_conditions[0]] > -1:
#         max_contr_2 = l_max_contr[condition] + r_max_contr[other_conditions[0]]
#         max_suppo_2 = l_max_suppo[condition] + r_max_suppo[other_conditions[0]]
#         max_oppos_2 = l_max_oppos[condition] + r_max_oppos[other_conditions[0]]
#
#     if l_max_contr[condition] > -1 and r_max_contr[other_conditions[1]] > -1:
#         max_contr_3 = l_max_contr[condition] + r_max_contr[other_conditions[1]]
#         max_suppo_3 = l_max_suppo[condition] + r_max_suppo[other_conditions[1]]
#         max_oppos_3 = l_max_oppos[condition] + r_max_oppos[other_conditions[1]]
#
#     if l_max_contr[condition] > -1 and r_max_contr[other_conditions[2]] > -1:
#         max_contr_4 = l_max_contr[condition] + r_max_contr[other_conditions[2]]
#         max_suppo_4 = l_max_suppo[condition] + r_max_suppo[other_conditions[2]]
#         max_oppos_4 = l_max_oppos[condition] + r_max_oppos[other_conditions[2]]
#
#     # special case: comparison with self
#     if l_max_contr[condition] > -1 and r_max_contr[condition] > -1:
#         max_contr_5 = l_max_contr[condition] + r_max_contr[condition]
#         max_suppo_5 = l_max_suppo[condition] + r_max_suppo[condition]
#         max_oppos_5 = l_max_oppos[condition] + r_max_oppos[condition]
#
#     # comparison with all others
#     if l_max_contr[4] > -1 and r_max_contr[condition] > -1:
#         max_contr_6 = l_max_contr[4] + r_max_contr[condition]
#         max_suppo_6 = l_max_suppo[4] + r_max_suppo[condition]
#         max_oppos_6 = l_max_oppos[4] + r_max_oppos[condition]
#
#     if l_max_contr[other_conditions[0]] > -1 and r_max_contr[condition] > -1:
#         max_contr_7 = l_max_contr[other_conditions[0]] + r_max_contr[condition]
#         max_suppo_7 = l_max_suppo[other_conditions[0]] + r_max_suppo[condition]
#         max_oppos_7 = l_max_oppos[other_conditions[0]] + r_max_oppos[condition]
#
#     if l_max_contr[other_conditions[1]] > -1 and r_max_contr[condition] > -1:
#         max_contr_8 = l_max_contr[other_conditions[1]] + r_max_contr[condition]
#         max_suppo_8 = l_max_suppo[other_conditions[1]] + r_max_suppo[condition]
#         max_oppos_8 = l_max_oppos[other_conditions[1]] + r_max_oppos[condition]
#
#     if l_max_contr[other_conditions[2]] > -1 and r_max_contr[condition] > -1:
#         max_contr_9 = l_max_contr[other_conditions[2]] + r_max_contr[condition]
#         max_suppo_9 = l_max_suppo[other_conditions[2]] + r_max_suppo[condition]
#         max_oppos_9 = l_max_oppos[other_conditions[2]] + r_max_oppos[condition]
#
#     max_contr = max(max_contr_1, max_contr_2, max_contr_3, max_contr_4, max_contr_5,
#                     max_contr_6, max_contr_7, max_contr_8, max_contr_9)
#
#     # Calculate maximum number of propairs, given a maxmimum number
#     # of pairs
#     max_suppo = -1
#     if max_contr == max_contr_1:
#         max_suppo = max(max_suppo, max_suppo_1)
#     if max_contr == max_contr_2:
#         max_suppo = max(max_suppo, max_suppo_2)
#     if max_contr == max_contr_3:
#         max_suppo = max(max_suppo, max_suppo_3)
#     if max_contr == max_contr_4:
#         max_suppo = max(max_suppo, max_suppo_4)
#     if max_contr == max_contr_5:
#         max_suppo = max(max_suppo, max_suppo_5)
#     if max_contr == max_contr_6:
#         max_suppo = max(max_suppo, max_suppo_6)
#     if max_contr == max_contr_7:
#         max_suppo = max(max_suppo, max_suppo_7)
#     if max_contr == max_contr_8:
#         max_suppo = max(max_suppo, max_suppo_8)
#     if max_contr == max_contr_9:
#         max_suppo = max(max_suppo, max_suppo_9)
#
#     # Calculate maximum number of antipairs, given a maxmimum number
#     # of pairs
#     max_oppos = -1
#     if max_contr == max_contr_1:
#         max_oppos = max(max_oppos, max_oppos_1)
#     if max_contr == max_contr_2:
#         max_oppos = max(max_oppos, max_oppos_2)
#     if max_contr == max_contr_3:
#         max_oppos = max(max_oppos, max_oppos_3)
#     if max_contr == max_contr_4:
#         max_oppos = max(max_oppos, max_oppos_4)
#     if max_contr == max_contr_5:
#         max_oppos = max(max_oppos, max_oppos_5)
#     if max_contr == max_contr_6:
#         max_oppos = max(max_oppos, max_oppos_6)
#     if max_contr == max_contr_7:
#         max_oppos = max(max_oppos, max_oppos_7)
#     if max_contr == max_contr_8:
#         max_oppos = max(max_oppos, max_oppos_8)
#     if max_contr == max_contr_9:
#         max_oppos = max(max_oppos, max_oppos_9)
#
#     return max_contr, max_suppo, max_oppos
#
#
# def _print_for_debugging(values):
#     pairs_max = pd.DataFrame(values[:, 0, :], columns=["11", "10", "01", "00", "0"])
#     pairs_supporting = pd.DataFrame(values[:, 1, :], columns=["11", "10", "01", "00", "0"])
#     pairs_opposing = pd.DataFrame(values[:, 2, :], columns=["11", "10", "01", "00", "0"])
# 
#     print('pairs_max', pairs_max, sep='\n')
#     print('pairs_supporting', pairs_supporting, sep='\n')
#     print('pairs_opposing', pairs_opposing, sep='\n')
