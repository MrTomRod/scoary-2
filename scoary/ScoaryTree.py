from __future__ import annotations

import logging
from typing import Union, Optional, Callable

import pandas as pd
from scipy.spatial import distance

from biotite.sequence.phylo import upgma
from biotite.sequence.phylo.tree import TreeNode as BiotiteTreeNode

from ete3 import Tree as EteTree

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

    def __init__(self, left: ScoaryTree = None, right: ScoaryTree = None, label: str = None):
        if left is None and right is None:
            self.is_leaf = True
            assert type(label) is str, f'A valid node has a label! {label=}'
            self.label = label
        else:
            self.is_leaf = False
            assert type(left) is ScoaryTree and type(right) is ScoaryTree, f'A valid tree has 0 or 2 children! {left=} {right=}'
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

    # def __copy__(self):
    #     if self.is_leaf:
    #         return ScoaryTree(label=self.label)
    #     else:
    #         return ScoaryTree(left=copy(self.left), right=copy(self.right))

    def rename(self, func: Callable):
        def convert(scoary_tree: ScoaryTree) -> ScoaryTree:
            """recursive function"""
            if scoary_tree.is_leaf:
                return ScoaryTree(label=func(scoary_tree.label))
            else:
                return ScoaryTree(left=convert(scoary_tree.left), right=convert(scoary_tree.right))

        return convert(self)

    @classmethod
    def from_newick(cls, newick: str) -> ScoaryTree:
        def convert(etetree: EteTree) -> ScoaryTree:
            """recursive function"""
            if len(etetree.children) == 0:
                return cls(label=etetree.name)
            else:
                return cls(left=convert(etetree.children[0]), right=convert(etetree.children[1]))

        ete_tree = EteTree(newick)
        return convert(ete_tree)

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
        def convert(node: BiotiteTreeNode) -> ScoaryTree:
            """recursive function"""
            if node.is_leaf():
                return cls(label=str(node_to_label[node]))
            else:
                return cls(left=convert(node.children[0]), right=convert(node.children[1]))

        distance_matrix = distance.squareform(distance.pdist(genes_df.T, 'hamming'))
        tree = upgma(distance_matrix)
        node_to_label: {BiotiteTreeNode: str} = {node: label for node, label in zip(tree.leaves, genes_df.columns)}
        return convert(tree.root)


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


def pick_contrasting(left_pairings: {(bool, bool): ScoaryTree}, right_pairings: {(bool, bool): ScoaryTree}) -> Optional[(ScoaryTree, ScoaryTree)]:
    return _pick([(True, True), (False, False), (True, False), (False, True)], left_pairings, right_pairings)


def pick_supporting(left_pairings, right_pairings) -> Optional[(ScoaryTree, ScoaryTree)]:
    return _pick([(True, True), (False, False)], left_pairings, right_pairings)


def pick_opposing(left_pairings, right_pairings) -> Optional[(ScoaryTree, ScoaryTree)]:
    return _pick([(True, False), (False, True)], left_pairings, right_pairings)


def count_max_pairings(tree: ScoaryTree, label_to_trait: {str: bool}, label_to_gene: {str: bool}, type: str) -> int:
    if type == 'contrasting':
        pick = pick_contrasting
    elif type == 'supporting':
        pick = pick_supporting
    elif type == 'opposing':
        pick = pick_opposing
    else:
        raise AssertionError(f"{type=} not in ('contrasting', 'supporting', 'opposing')")

    n_picks = 0

    logger.info(f'## pick {type} pairs')

    def pick_pairings(scoary_tree: ScoaryTree) -> {(bool, bool): ScoaryTree}:
        if scoary_tree.is_leaf:
            return {(label_to_gene[scoary_tree.label], label_to_trait[scoary_tree.label]): scoary_tree}
        else:
            left_pairings = pick_pairings(scoary_tree.left)
            right_pairings = pick_pairings(scoary_tree.right)
            pair = pick(left_pairings, right_pairings)
            if pair is not None:
                logger.info(f"{type} pick: {pair}")
                nonlocal n_picks
                n_picks += 1
                return {}
            else:
                return left_pairings | right_pairings

    pick_pairings(tree)

    return n_picks
