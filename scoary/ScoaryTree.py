from __future__ import annotations

from functools import cached_property
from typing import Optional, Callable

import numpy as np
import pandas as pd
from scipy.spatial import distance

from .newick import parse_newick
from .upgma import upgma


class ScoaryTree:
    left: Optional[ScoaryTree] = None
    right: Optional[ScoaryTree] = None
    label: Optional[str] = None
    is_leaf: bool = False
    _values: Optional[np.ndarray] = None
    _prune = False

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
        return self._newick()

    def __repr__(self):
        return str(self)

    def _newick(self):
        return self.label if self.is_leaf else f"({self.left._newick()},{self.right._newick()})"

    def to_newick(self) -> str:
        return f'{self._newick()};'

    def write_newick(self, path: str):
        with open(path, 'w') as f:
            f.write(self.to_newick())

    def labels(self) -> [str]:
        if self.is_leaf:
            return [self.label]
        else:
            return self.left.labels() + self.right.labels()

    def uniquify(self, label_to_trait: {str: bool}):
        def uniquify(tree: ScoaryTree) -> str:
            if tree.is_leaf:
                return '1' if label_to_trait[tree.label] else '0'
            else:
                l, r = uniquify(tree.left), uniquify(tree.right)
                return f'({l}{r})' if l < r else f'({r}{l})'

        return uniquify(self)

    def copy(self):
        def copy(tree: ScoaryTree) -> ScoaryTree:
            if tree.is_leaf:
                return ScoaryTree(label=tree.label)
            else:
                return ScoaryTree(left=copy(tree.left), right=copy(tree.right))

        return copy(self)

    def prune(self, labels: [str]) -> ScoaryTree:
        n_labels_found = 0

        def prune(tree: ScoaryTree) -> Optional[ScoaryTree]:
            if tree.is_leaf:
                if tree.label in labels:
                    nonlocal n_labels_found
                    n_labels_found += 1
                    return ScoaryTree(label=tree.label)
                else:
                    return None
            else:
                left, right = prune(tree.left), prune(tree.right)
                if left and right:
                    return ScoaryTree(left=left, right=right)
                if left:
                    return left
                if right:
                    return right
                return None

        pruned_tree = prune(self)

        if n_labels_found != len(labels):
            missing = set(labels).difference(set(self.labels()))
            raise AssertionError(f'Pruning went wrong: did not find all labels in tree! '
                                 f'{n_labels_found=}; {missing=}; tree={self}')

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
                    this._prune = True  # mark for pruning

                # append self to parent
                setattr(current_parent, current_direction, this)

                if current_direction == 'right':
                    # prune
                    current_parent.__prune()
                    stack.pop()

                    # found terminal node
                    # # GO UP UNTIL CAN GO RIGHT
                    while stack:
                        ancestor_node, ancestor_direction = stack[-1]
                        if ancestor_direction == 'right':
                            ancestor_node.__prune()
                            stack.pop()
                        else:
                            break

                    if not stack:
                        print(f'done\n{self}\n{root}')
                        break

                # pop left node -> go right next
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

                if current_direction == 'right':
                    # found terminal node
                    # # GO UP UNTIL CAN GO RIGHT
                    while stack and stack[-1][1] == 'right':
                        stack.pop()

                    if not stack:
                        print(f'done\n{self}\n{root}')
                        break

                # pop left node -> go right next
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

        def convert(tree: ScoaryTree) -> ScoaryTree:
            """recursive function"""
            if tree.is_leaf:
                return ScoaryTree(label=func(tree.label))
            else:
                return ScoaryTree(left=convert(tree.left), right=convert(tree.right))

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

    @cached_property
    def to_list(self) -> []:
        def to_list(tree: ScoaryTree) -> str | []:
            if tree.is_leaf:
                return tree.label
            else:
                return [to_list(tree.left), to_list(tree.right)]

        return to_list(self)

    @classmethod
    def from_presence_absence(cls, genes_df: pd.DataFrame) -> ScoaryTree:
        distance_matrix = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')),
                                       columns=genes_df.columns)
        tree_as_list = upgma(distance_matrix)
        return cls.from_list(tree_as_list)

    def __prune(self):
        if self.left._prune and self.right._prune:
            self._prune = True
        elif self.left._prune:
            # become right
            self.label = self.right.label
            self.is_leaf = self.right.is_leaf
            self.left = self.right.left
            self.right = self.right.right
        elif self.right._prune:
            # become left
            self.label = self.left.label
            self.is_leaf = self.left.is_leaf
            self.right = self.left.right
            self.left = self.left.left
