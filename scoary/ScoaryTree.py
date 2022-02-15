from __future__ import annotations

import logging
from typing import Optional, Callable

import pandas as pd
from scipy.spatial import distance

from .newick import parse_newick
from .upgma import upgma

# from numba import njit, prange

logger = logging.getLogger('scoary-tree')


def union_dict_of_lists(*dicts: [{str: list}]) -> {str: list}:
    """
    Combines multiple dictionaries where the value is a list by combining the lists.

    :param dicts: list of dictionaries that have lists as values
    :return: combined dictionary
    """
    res = {}
    for dict in dicts:
        for key, list_ in dict.items():
            if key in res:
                res[key].extend(list_)
            else:
                res[key] = list_
    return res


class ScoaryTree:
    left: Optional[ScoaryTree] = None
    right: Optional[ScoaryTree] = None
    label: Optional[str] = None
    is_leaf: bool = False
    __prune = False

    def __init__(self, left: ScoaryTree = None, right: ScoaryTree = None, label: str = None):
        if left is None and right is None:
            self.is_leaf = True
            assert type(label) is str, f'A valid node has a label! {label=}'
            self.label = label
        else:
            self.is_leaf = False
            assert type(left) is ScoaryTree and type(
                right) is ScoaryTree, f'A valid tree has 0 or 2 children! {left=} {right=}'
            self.left = left
            self.right = right

    def __str__(self) -> str:
        return self.label if self.is_leaf else f"({self.left},{self.right})"

    def __repr__(self):
        return str(self)

    def to_newick(self) -> str:
        return f'{self};'

    def labels(self) -> [str]:
        if self.is_leaf:
            return [self.label]
        else:
            return self.left.labels() + self.right.labels()

    def copy(self):
        def copy(scoary_tree: ScoaryTree) -> ScoaryTree:
            if scoary_tree.is_leaf:
                return ScoaryTree(label=scoary_tree.label)
            else:
                return ScoaryTree(left=copy(scoary_tree.left), right=copy(scoary_tree.right))

        return copy(self)

    def prune(self, labels: [str]) -> ScoaryTree:
        n_labels_found = 0

        def prune(scoary_tree: ScoaryTree) -> Optional[ScoaryTree]:
            if scoary_tree.is_leaf:
                if scoary_tree.label in labels:
                    nonlocal n_labels_found
                    n_labels_found += 1
                    return ScoaryTree(label=scoary_tree.label)
                else:
                    return None
            else:
                left, right = prune(scoary_tree.left), prune(scoary_tree.right)
                if left and right:
                    return ScoaryTree(left=left, right=right)
                if left:
                    return left
                if right:
                    return right
                return None

        pruned_tree = prune(self)

        assert n_labels_found == len(
            labels), f'Pruning went wrong: did not find all labels in tree! {n_labels_found=}; {labels=}; tree={self}'

        return pruned_tree

    def prune_nonrecursive(self, labels: [str]) -> ScoaryTree:
        if self.is_leaf:
            assert [self.label] == labels, f'Pruning went wrong. {[self.label]} != {labels}'
            return ScoaryTree(label=self.label)

        n_labels_found = 0

        root = ScoaryTree(left=self.left, right=self.right)

        stack = [(root, 'right'), (root, 'left')]

        while stack:
            current_parent, current_direction = stack[-1]
            current_node: ScoaryTree = getattr(current_parent, current_direction)

            if current_node.is_leaf:
                # current node is leaf

                this = ScoaryTree(label=current_node.label)

                if current_node.label in labels:
                    n_labels_found += 1
                else:
                    this.__prune = True  # mark for pruning

                # append self to parent
                setattr(current_parent, current_direction, this)

                if current_direction == 'left':
                    # remove leaf from stack
                    stack.pop()
                else:
                    # prune
                    current_parent._prune()
                    stack.pop()

                    # found terminal node
                    # # GO UP UNTIL CAN GO RIGHT
                    while stack:
                        ancestor_node, ancestor_direction = stack[-1]
                        if ancestor_direction == 'right':
                            ancestor_node._prune()
                            stack.pop()
                        else:
                            break

                    if not stack:
                        print(f'done\n{self}\n{root}')
                        break

                    # pop one more -> go right on this node
                    stack.pop()

            else:
                this = ScoaryTree(left=current_node.left, right=current_node.right)
                stack.extend([(this, 'right'), (this, 'left')])

                # append self to parent
                setattr(current_parent, current_direction, this)

        return root

    def copy_nonrecursive(self) -> ScoaryTree:
        return self.rename_nonrecursive(func=lambda label: label)

    def rename_nonrecursive(self, func: Callable):
        if self.is_leaf:
            return ScoaryTree(label=func(self.label))

        root = ScoaryTree(left=self.left, right=self.right)

        stack = [(root, 'right'), (root, 'left')]

        while stack:
            current_parent, current_direction = stack[-1]
            current_node: ScoaryTree = getattr(current_parent, current_direction)

            if current_node.is_leaf:
                # current node is leaf
                this = ScoaryTree(label=func(current_node.label))

                # append self to parent
                setattr(current_parent, current_direction, this)

                # remove leaf from stack
                stack.pop()

                if current_direction == 'right':
                    # found terminal node
                    # # GO UP UNTIL CAN GO RIGHT
                    while stack and stack[-1][1] == 'right':
                        stack.pop()

                    if not stack:
                        print(f'done\n{self}\n{root}')
                        break

                    # pop one more -> go right on this node
                    stack.pop()

            else:
                this = ScoaryTree(left=current_node.left, right=current_node.right)
                stack.extend([(this, 'right'), (this, 'left')])

                # append self to parent
                setattr(current_parent, current_direction, this)

        return root

    def rename(self, func: Callable):
        """
        Apply a function to each leaf label.

        Only used for debugging. This recursive function could cause RecursionError for big trees.
        """

        def convert(scoary_tree: ScoaryTree) -> ScoaryTree:
            """recursive function"""
            if scoary_tree.is_leaf:
                return ScoaryTree(label=func(scoary_tree.label))
            else:
                return ScoaryTree(left=convert(scoary_tree.left), right=convert(scoary_tree.right))

        return convert(self)

    @classmethod
    def from_newick(cls, newick: str) -> ScoaryTree:
        list_tree = parse_newick(newick)
        return cls.from_list(list_tree)

    @classmethod
    def from_list(cls, tree: []) -> ScoaryTree:
        def convert(list_tree):
            """recursive function"""
            if type(list_tree) is str:
                return cls(label=list_tree)
            else:
                return cls(left=convert(list_tree[0]), right=convert(list_tree[1]))

        return convert(tree)

    def to_list(self) -> []:
        if self.is_leaf:
            return self.label
        else:
            return [self.left.to_list(), self.right.to_list()]

    @classmethod
    def from_presence_absence(cls, genes_df: pd.DataFrame) -> ScoaryTree:
        distance_matrix = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')),
                                       columns=genes_df.columns)
        tree_as_list = upgma(distance_matrix)
        return cls.from_list(tree_as_list)

    def _prune(self):
        if self.left.__prune and self.right.__prune:
            self.__prune = True
        elif self.left.__prune:
            # become right
            self.label = self.right.label
            self.is_leaf = self.right.is_leaf
            self.left = self.right.left
            self.right = self.right.right
        elif self.right.__prune:
            # become left
            self.label = self.left.label
            self.is_leaf = self.left.is_leaf
            self.right = self.left.right
            self.left = self.left.left


def _pick(
        pairs_to_pick: [(bool, bool)],
        left_pairings: {(bool, bool): ScoaryTree},
        right_pairings: {(bool, bool): ScoaryTree}
) -> Optional[(ScoaryTree, ScoaryTree)]:
    for pair in pairs_to_pick:
        antipair = (not pair[0], not pair[1])
        if pair in left_pairings and antipair in right_pairings:
            return left_pairings[pair], right_pairings[antipair]
    return None  # found no such pairing

# def pick_contrasting(left_pairings: {(bool, bool): ScoaryTree}, right_pairings: {(bool, bool): ScoaryTree}) -> Optional[(ScoaryTree, ScoaryTree)]:
#     return _pick([(True, True), (False, False), (True, False), (False, True)], left_pairings, right_pairings)
#
#
# def pick_supporting(left_pairings, right_pairings) -> Optional[(ScoaryTree, ScoaryTree)]:
#     return _pick([(True, True), (False, False)], left_pairings, right_pairings)
#
#
# def pick_opposing(left_pairings, right_pairings) -> Optional[(ScoaryTree, ScoaryTree)]:
#     return _pick([(True, False), (False, True)], left_pairings, right_pairings)
#
#
# def count_max_pairings(tree: ScoaryTree, label_to_trait: {str: bool}, label_to_gene: {str: bool}, type: str) -> int:
#     if type == 'contrasting':
#         pick = pick_contrasting
#     elif type == 'supporting':
#         pick = pick_supporting
#     elif type == 'opposing':
#         pick = pick_opposing
#     else:
#         raise AssertionError(f"{type=} not in ('contrasting', 'supporting', 'opposing')")
#
#     n_picks = 0
#
#     logger.info(f'## pick {type} pairs')
#
#     def pick_pairings(scoary_tree: ScoaryTree) -> {(bool, bool): ScoaryTree}:
#         if scoary_tree.is_leaf:
#             return {(label_to_gene[scoary_tree.label], label_to_trait[scoary_tree.label]): scoary_tree}
#         else:
#             left_pairings = pick_pairings(scoary_tree.left)
#             right_pairings = pick_pairings(scoary_tree.right)
#             pair = pick(left_pairings, right_pairings)
#             if pair is not None:
#                 logger.info(f"{type} pick: {pair}")
#                 nonlocal n_picks
#                 n_picks += 1
#                 return {}
#             else:
#                 return left_pairings | right_pairings
#
#     pick_pairings(tree)
#
#     return n_picks
#
#
# def count_sup_op(tree: ScoaryTree, label_to_trait: {str: bool}, label_to_gene: {str: bool}) -> (int, int):
#     n_sup = 0
#     n_opp = 0
#
#     logger.info(f'## pick count_sup_op pairs')
#
#     def pick_pairings(scoary_tree: ScoaryTree) -> {(bool, bool): ScoaryTree}:
#         if scoary_tree.is_leaf:
#             return {(label_to_gene[scoary_tree.label], label_to_trait[scoary_tree.label]): scoary_tree}
#         else:
#             left_pairings = pick_pairings(scoary_tree.left)
#             right_pairings = pick_pairings(scoary_tree.right)
#             supp_pair = pick_supporting(left_pairings, right_pairings)
#             if supp_pair is not None:
#                 logger.info(f"supp pick: {supp_pair}")
#                 nonlocal n_sup
#                 n_sup += 1
#                 return {}
#             else:
#                 opp_pair = pick_opposing(left_pairings, right_pairings)
#                 if opp_pair is not None:
#                     logger.info(f"opp pick: {opp_pair}")
#                     nonlocal n_opp
#                     n_opp += 1
#                     return {}
#                 else:
#                     return left_pairings | right_pairings
#
#     pick_pairings(tree)
#
#     return n_sup, n_opp
#
#
# def count_best_worst(tree: ScoaryTree, label_to_trait: {str: bool}, label_to_gene: {str: bool}) -> (int, int):
#     n_best = 0
#     n_worst = 0
#
#     logger.info(f'## pick count_sup_op pairs')
#
#     def pick_pairings(scoary_tree: ScoaryTree) -> {(bool, bool): ScoaryTree}:
#         if scoary_tree.is_leaf:
#             return {(label_to_gene[scoary_tree.label], label_to_trait[scoary_tree.label]): scoary_tree}
#         else:
#             left_pairings = pick_pairings(scoary_tree.left)
#             right_pairings = pick_pairings(scoary_tree.right)
#
#             supp_pair = pick_supporting(left_pairings, right_pairings)
#             opp_pair = pick_opposing(left_pairings, right_pairings)
#
#             if supp_pair is None and opp_pair is None:
#                 return left_pairings | right_pairings
#
#             if supp_pair is not None:
#                 logger.info(f"supp pick: {supp_pair}")
#                 nonlocal n_best
#                 n_best += 1
#
#             if opp_pair is not None:
#                 logger.info(f"opp pick: {opp_pair}")
#                 nonlocal n_worst
#                 n_worst += 1
#
#             return {}
#
#     pick_pairings(tree)
#
#     return n_best, n_worst
