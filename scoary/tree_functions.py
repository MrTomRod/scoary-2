from typing import Union

import pandas as pd
from scipy.spatial import distance

from biotite.sequence.phylo import upgma
from biotite.sequence.phylo.tree import TreeNode, Tree

from ete3 import Tree as EteTree


def _biotite_to_list(tree: Tree, labels: [str]) -> Union[str, list]:
    node_to_label: {TreeNode: str} = {node: label for node, label in zip(tree.leaves, labels)}

    def _tree_to_list(node: TreeNode):
        if node.is_leaf():
            return str(node_to_label[node])
        else:
            return [_tree_to_list(child) for child in node.children]

    return _tree_to_list(tree.root)


def _ete_to_list(etetree: EteTree) -> Union[str, list]:
    if len(etetree.children) == 0:
        return etetree.name
    else:
        return [_ete_to_list(node) for node in etetree.children]


def tree_from_presence_absence(genes_df: pd.DataFrame) -> []:
    distance_matrix = distance.squareform(distance.pdist(genes_df.T, 'hamming'))
    tree = upgma(distance_matrix)
    return _biotite_to_list(tree, labels=genes_df.columns)


def tree_from_file(newick: str) -> []:
    return _ete_to_list(EteTree(newick))


def tree_to_newick(tree: []) -> str:
    def newickfiy(tree):
        if type(tree) is str:
            return tree
        else:
            return f"({','.join(newickfiy(t) for t in tree)})"

    return f'{newickfiy(tree)};'
