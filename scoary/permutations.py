import os
import logging
import pandas as pd
import numpy as np
import scipy.stats as ss

from .KeyValueStore import KeyValueStore
from .picking import pick
from .ScoaryTree import ScoaryTree

logger = logging.getLogger('scoary.permutations')


class ConfintStore(KeyValueStore):
    def create_db(self):
        self._create_db(
            columns={
                'tree': 'str',
                'n_pos_assoc': 'int',
                'n_permut': 'int',
                'confidence_interval': 'str'
            },
            pk_col='tree, n_pos_assoc, n_permut'
        )

    def get(self, tree: str, n_pos_assoc: int, n_permut: int):
        sql = f'SELECT confidence_interval FROM {self.table_name} WHERE tree = ? AND n_pos_assoc = ? AND n_permut = ?'
        res = self.cur.execute(
            sql,
            (tree, n_pos_assoc, n_permut,)
        ).fetchone()
        return np.frombuffer(res[0], dtype=float) if res is not None else None

    def set(self, tree: str, n_pos_assoc: int, n_permut: int, confidence_interval: [float]):
        confidence_interval = confidence_interval.tobytes()
        sql = f'INSERT OR IGNORE INTO {self.table_name} VALUES (?, ?, ?, ?)'
        self.cur.execute(
            sql, (tree, n_pos_assoc, n_permut, confidence_interval)
        )
        self.con.commit()


CONFINT_CACHE = ConfintStore(table_name='confint_cache', db_path=os.environ.get('CONFINT_DB', None))


def create_permuted_df(labels: [str], n_positive: int, n_permut: int, random_state: int = None):
    if random_state:
        np.random.seed(random_state)

    n_negative = len(labels) - n_positive
    arr = np.repeat(np.array([[1] * n_positive + [0] * n_negative]), n_permut, axis=0)

    # creates a copy -> slow
    arr = np.apply_along_axis(np.random.permutation, axis=1, arr=arr)

    return pd.DataFrame(arr, columns=labels)


def permute_picking(
        trait: str,
        tree: ScoaryTree,
        label_to_trait: pd.Series | dict,
        result_df: pd.DataFrame,
        genes_bool_df: pd.DataFrame,
        n_permut: int,
        random_state: int = None,
) -> np.array:
    if type(label_to_trait) is dict:
        label_to_trait = pd.Series(label_to_trait, dtype='boolean')
    n_tot = len(label_to_trait)
    n_pos = sum(label_to_trait)
    n_neg = n_tot - n_pos
    labels = label_to_trait.keys()

    n_reused = 0

    pvals = []
    for _, row in result_df.iterrows():
        label_to_gene = genes_bool_df.loc[row.Gene]
        unique_topology = tree.uniquify(label_to_gene)

        is_positively_correlated = row.supporting >= row.opposing
        estimator = (row.supporting if is_positively_correlated else row.opposing) / row.contrasting
        n_pos_assoc = n_pos if is_positively_correlated else n_neg

        permuted_estimators = CONFINT_CACHE.get(unique_topology, n_pos_assoc, n_permut)
        if permuted_estimators is None:
            permuted_df = create_permuted_df(
                labels=labels, n_positive=n_pos_assoc,
                n_permut=n_permut, random_state=random_state
            )
            max_contr, max_suppo, max_oppos = pick(
                tree=tree.to_list, label_to_trait_a=label_to_gene,
                trait_b_df=permuted_df, calc_pvals=False
            )

            permuted_estimators = max_suppo / max_contr
            CONFINT_CACHE.set(unique_topology, n_pos_assoc, n_permut, permuted_estimators)
        else:
            n_reused += 1

        pval = ((permuted_estimators >= estimator).sum() + 1) / (n_permut + 1)
        pvals.append(pval)

    logger.debug(f'{trait}: reused {n_reused} out of {len(result_df)}')

    return pvals
