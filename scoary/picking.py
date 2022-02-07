import os

from numba import jit, njit, prange, boolean

import numpy as np
import pandas as pd

PRINT_ANY = os.environ.get('PRINT_ANY', 'FALSE').lower() == 'true'
PRINT_INIT = os.environ.get('PRINT_INIT', 'FALSE').lower() == 'true'
PRINT_COMBINE = os.environ.get('PRINT_COMBINE', 'FALSE').lower() == 'true'
PRINT_BEFORE_COMBINE = os.environ.get('PRINT_BEFORE_COMBINE', 'FALSE').lower() == 'true'
PRINT_AFTER_COMBINE = os.environ.get('PRINT_AFTER_COMBINE', 'FALSE').lower() == 'true'


def pick(tree: [], label_to_gene: {str: bool}, permuted_traits_df: pd.DataFrame):
    def _pick(left_label, right_label):
        if type(left_label) is str:
            if PRINT_ANY: print(f'INIT LEAF LEFT: {left_label}')
            left = init_leaf(
                gene=label_to_gene[left_label],
                traits=permuted_traits_df[left_label].to_numpy()
            )
        else:
            left = _pick(left_label[0], left_label[1])
        if type(right_label) is str:
            if PRINT_ANY: print(f'INIT LEAF RIGHT: {right_label}')
            right = init_leaf(
                gene=label_to_gene[right_label],
                traits=permuted_traits_df[right_label].to_numpy()
            )
        else:
            right = _pick(right_label[0], right_label[1])

        if PRINT_ANY: print(f'COMBINE BRANCHES [{left_label}, {right_label}]')
        combined = combine_branches(left, right)
        if PRINT_COMBINE:
            print(combined)
        return combined

    values = _pick(tree[0], tree[1])

    max_contrasting = values[:, 0, :].max(axis=1)
    max_supporting = values[:, 1, :].max(axis=1)
    max_opposing = values[:, 2, :].max(axis=1)

    assert len(max_contrasting) == len(max_supporting) == len(max_opposing) == len(permuted_traits_df)

    return max_contrasting, max_supporting, max_opposing


# selecting:values[<TRAITS>, <3 TYPES OF PAIRINGS>, <5 COMBINATIONS>]
# selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: AB; 1: Ab; 2: aB; 3: ab; 4: free>]

# all pairs_max: values[:, 0, :]
# all pairs_supporting: values[:, 1, :]
# all pairs_opposing: values[:, 2, :]

# selecting first trait: values[0,:,:]
# selecting all AB: values[:, 0, :]

# @njit('int64[:, ::3, ::5](b1, boolean[:])', nogil=True, boundscheck=False)
def init_leaf(gene: bool, traits: np.array) -> np.array:
    n_traits = traits.shape[0]

    values = np.full(shape=(n_traits, 3, 5), fill_value=-1, dtype='int')
    if gene:
        for i, trait in enumerate(traits):
            if trait:
                values[i, 0, 0] = 0  # pairs_max
                values[i, 1, 0] = 0  # pairs_supporting
                values[i, 2, 0] = 0  # pairs_opposing
            else:
                values[i, 0, 1] = 0  # pairs_max
                values[i, 1, 1] = 0  # pairs_supporting
                values[i, 2, 1] = 0  # pairs_opposing
    else:
        for i, trait in enumerate(traits):
            if trait:
                values[i, 0, 2] = 0  # pairs_max
                values[i, 1, 2] = 0  # pairs_supporting
                values[i, 2, 2] = 0  # pairs_opposing
            else:
                values[i, 0, 3] = 0  # pairs_max
                values[i, 1, 3] = 0  # pairs_supporting
                values[i, 2, 3] = 0  # pairs_opposing

    # print_for_debugging(values)
    if PRINT_INIT:
        print(values)

    return values


def print_for_debugging(values):
    pairs_max = pd.DataFrame(values[:, 0, :], columns=["AB", "Ab", "aB", "ab", "0"])
    pairs_supporting = pd.DataFrame(values[:, 1, :], columns=["AB", "Ab", "aB", "ab", "0"])
    pairs_opposing = pd.DataFrame(values[:, 2, :], columns=["AB", "Ab", "aB", "ab", "0"])

    print('pairs_max', pairs_max, sep='\n')
    print('pairs_supporting', pairs_supporting, sep='\n')
    print('pairs_opposing', pairs_opposing, sep='\n')


# @njit('int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False)
def calculate_max_nofree(left: np.array, right: np.array):
    values = np.full(shape=(3, 5), fill_value=-1, dtype='int')

    if left[0][4] > -1 and right[0][4] > -1:  # 0 vs 0
        values[0][0] = left[0][4] + right[0][4]
        values[1][0] = left[1][4] + right[1][4]
        values[2][0] = left[2][4] + right[2][4]

    if left[0][0] > -1 and right[0][3] > -1:  # AB vs ab
        values[0][1] = left[0][0] + right[0][3] + 1
        values[1][1] = left[1][0] + right[1][3] + 1
        values[2][1] = left[2][0] + right[2][3]

    if left[0][3] > -1 and right[0][0] > -1:  # ab vs AB
        values[0][2] = left[0][3] + right[0][0] + 1
        values[1][2] = left[1][3] + right[1][0] + 1
        values[2][2] = left[2][3] + right[2][0]

    if left[0][1] > -1 and right[0][2] > -1:  # Ab vs aB
        values[0][3] = left[0][1] + right[0][2] + 1
        values[1][3] = left[1][1] + right[1][2]
        values[2][3] = left[2][1] + right[2][2] + 1

    if left[0][2] > -1 and right[0][1] > -1:  # aB vs Ab
        values[0][4] = left[0][2] + right[0][1] + 1
        values[1][4] = left[1][2] + right[1][1]
        values[2][4] = left[2][2] + right[2][1] + 1

    max_contrasting = values[0].max()

    max_supporting = -1
    for i in range(5):
        if values[0][i] == max_contrasting and values[1][i] > max_supporting:
            max_supporting = values[1][i]

    max_opposing = -1
    for i in range(5):
        if values[0][i] == max_contrasting and values[2][i] > max_opposing:
            max_opposing = values[2][i]

    return max_contrasting, max_supporting, max_opposing


# @njit('int64, int64[::3, ::5], int64[::3, ::5]', nogil=True, boundscheck=False)
def calculate_max_condition(condition, left: np.array, right: np.array):
    values = np.full(shape=(3, 9), fill_value=-1, dtype='int')

    other_conditions = [i for i in range(4) if i != condition]

    if left[0][condition] > -1:
        # comparison with other_conditions
        if right[0][other_conditions[0]] > -1:
            values[0][0] = left[0][condition] + right[0][other_conditions[0]]
            values[1][0] = left[1][condition] + right[1][other_conditions[0]]
            values[2][0] = left[2][condition] + right[2][other_conditions[0]]

        if right[0][other_conditions[1]] > -1:
            values[0][1] = left[0][condition] + right[0][other_conditions[1]]
            values[1][1] = left[1][condition] + right[1][other_conditions[1]]
            values[2][1] = left[2][condition] + right[2][other_conditions[1]]

        if right[0][other_conditions[2]] > -1:
            values[0][2] = left[0][condition] + right[0][other_conditions[2]]
            values[1][2] = left[1][condition] + right[1][other_conditions[2]]
            values[2][2] = left[2][condition] + right[2][other_conditions[2]]

        # special case: comparison with self
        if right[0][condition] > -1:
            values[0][3] = left[0][condition] + right[0][condition]
            values[1][3] = left[1][condition] + right[1][condition]
            values[2][3] = left[2][condition] + right[2][condition]

        if right[0][4] > -1:
            values[0][4] = left[0][condition] + right[0][4]
            values[1][4] = left[1][condition] + right[1][4]
            values[2][4] = left[2][condition] + right[2][4]


    if right[0][condition] > -1:
        # comparison with all others
        if left[0][4] > -1:
            values[0][5] = left[0][4] + right[0][condition]
            values[1][5] = left[1][4] + right[1][condition]
            values[2][5] = left[2][4] + right[2][condition]

        if left[0][other_conditions[0]] > -1:
            values[0][6] = left[0][other_conditions[0]] + right[0][condition]
            values[1][6] = left[1][other_conditions[0]] + right[1][condition]
            values[2][6] = left[2][other_conditions[0]] + right[2][condition]

        if left[0][other_conditions[1]] > -1:
            values[0][7] = left[0][other_conditions[1]] + right[0][condition]
            values[1][7] = left[1][other_conditions[1]] + right[1][condition]
            values[2][7] = left[2][other_conditions[1]] + right[2][condition]

        if left[0][other_conditions[2]] > -1:
            values[0][8] = left[0][other_conditions[2]] + right[0][condition]
            values[1][8] = left[1][other_conditions[2]] + right[1][condition]
            values[2][8] = left[2][other_conditions[2]] + right[2][condition]

    max_contrasting = values[0].max()

    max_supporting = -1
    for i in range(9):
        if values[0][i] == max_contrasting and values[1][i] > max_supporting:
            max_supporting = values[1][i]

    max_opposing = -1
    for i in range(9):
        if values[0][i] == max_contrasting and values[2][i] > max_opposing:
            max_opposing = values[2][i]

    return max_contrasting, max_supporting, max_opposing


# @njit('int64[:, ::3, ::5], int64[:, ::3, ::5]', nogil=True, boundscheck=False)  # , parallel=True)
def combine_branches(left: np.array, right: np.array):
    assert left.shape == right.shape
    n_traits = left.shape[0]

    if PRINT_BEFORE_COMBINE:
        print('BEFORE_COMBINE_L\n', left)
        print('BEFORE_COMBINE_R\n', right)

    values = np.full(shape=left.shape, fill_value=-1, dtype='int')

    # selecting:values[<TRAITS>, <0: max; 1: supporting; 2: opposing>, <0: AB; 1: Ab; 2: aB; 3: ab; 4: free>]
    for trait_id in prange(n_traits):
        for cond in prange(4):  # {"AB": 0, "Ab": 1, "aB": 2, "ab": 3, "0": 4}
            max_contrasting, max_supporting, max_opposing = calculate_max_condition(
                cond,
                left[trait_id, :, :],
                right[trait_id, :, :]
            )
            values[trait_id, 0, cond] = max_contrasting
            values[trait_id, 1, cond] = max_supporting
            values[trait_id, 2, cond] = max_opposing
        max_contrasting, max_supporting, max_opposing = calculate_max_nofree(left[trait_id, :, :],
                                                                             right[trait_id, :, :])
        values[trait_id, 0, 4] = max_contrasting
        values[trait_id, 1, 4] = max_supporting
        values[trait_id, 2, 4] = max_opposing

        if PRINT_AFTER_COMBINE:
            print('AFTER_COMBINE\n', values)
    return values
