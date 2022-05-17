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
        trait:str,
        tree: ScoaryTree,
        label_to_trait: {str: bool},
        result_df: pd.DataFrame,
        all_label_to_gene: dict[str:bool],
        n_permut: int,
        random_state: int = None
) -> np.array:
    n_total = len(label_to_trait)
    n_trait = sum(label_to_trait.values())
    n_no_trait = n_total - n_trait
    labels = label_to_trait.keys()

    n_reused = 0

    pvals = []
    for _, row in result_df.iterrows():
        label_to_gene = all_label_to_gene[row.Gene]
        unique_topology = tree.uniquify(label_to_gene)

        is_positively_correlated = row.supporting >= row.opposing
        estimator = (row.supporting if is_positively_correlated else row.opposing) / row.contrasting
        n_pos_assoc = n_trait if is_positively_correlated else n_no_trait

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


def permute_gene_picking(
        tree: ScoaryTree,
        label_to_trait: {str: bool},
        result_df: pd.DataFrame,
        all_label_to_gene: dict[str:bool],
        n_permut: int
) -> np.array:
    unique_string = tree.uniquify(label_to_trait)
    if unique_string not in CONFINT_CACHE:
        CONFINT_CACHE[unique_string] = {}
    CURRENT_CACHE = CONFINT_CACHE[unique_string]

    n_total = len(label_to_trait)
    labels = label_to_trait.keys()

    n_reused = 0

    pvals = []
    for _, row in result_df.iterrows():
        label_to_gene = all_label_to_gene[row.Gene]
        n_genes = sum(label_to_gene.values())

        is_positively_correlated = row.supporting >= row.opposing
        estimator = (row.supporting if is_positively_correlated else row.opposing) / row.contrasting
        n_pos_assoc = n_genes if is_positively_correlated else n_total - n_genes

        if n_pos_assoc not in CURRENT_CACHE:
            permuted_df = create_permuted_df(labels=labels, n_positive=n_pos_assoc, n_permut=n_permut)
            max_contr, max_suppo, max_oppos = pick(
                tree=tree.to_list, label_to_trait_a=label_to_trait,
                trait_b_df=permuted_df, calc_pvals=False
            )

            permuted_estimators = max_suppo / max_contr
            CURRENT_CACHE[n_pos_assoc] = permuted_estimators
        else:
            n_reused += 1
            logger.debug(f'reusing {row.Gene}')
            permuted_estimators = CURRENT_CACHE[n_pos_assoc]

        pval = ((permuted_estimators >= estimator).sum() + 1) / (n_permut + 1)
        pvals.append(pval)

    logger.debug(f'reused {n_reused} out of {len(result_df)}')

    return pvals


def permute_trait_picking_scoary1(
        genes: pd.Series,
        all_label_to_gene: {str: bool},
        tree: ScoaryTree,
        label_to_trait: {str: bool},
        n_permut: int
) -> [int]:
    from .scoary_1_picking import convert_upgma_to_phylotree, permute_gtc
    tree = tree.prune(labels=label_to_trait.keys())
    labels = tree.labels()
    tree = tree.to_list
    boolify = lambda t1, t2: f"{'A' if t1 else 'a'}{'B' if t2 else 'b'}"
    pvals = []
    for gene in genes:
        label_to_gene = all_label_to_gene[gene]
        gtc = {l: boolify(label_to_gene[l], label_to_trait[l]) for l in labels}
        r = 0
        _MyPhyloTree_, Unpermuted_tree = convert_upgma_to_phylotree(tree, gtc)
        Pro = 'Pro' if Unpermuted_tree['Pro'] >= Unpermuted_tree['Anti'] else 'Anti'
        Unpermuted_estimator = (float(Unpermuted_tree[Pro]) / Unpermuted_tree["Total"])
        for i in range(n_permut):
            PermutedGTC = permute_gtc(gtc)
            _MyPhyloTree_, NewPhyloTree = convert_upgma_to_phylotree(tree, PermutedGTC)
            if (float(NewPhyloTree[Pro]) / NewPhyloTree["Total"] >=
                    Unpermuted_estimator):
                r += 1

            # Check how many estimators are higher than the unpermuted
            # If, after more than 30 iterations, r indicates a p > 0.1,
            # abort
            if i >= 30:
                if (1 - ss.binom.cdf(r, i, 0.1)) < 0.05:
                    emp_p = (r + 1.0) / (i + 2.0)
                    print(f'break {gene} {emp_p}')
                    break
        else:
            emp_p = (r + 1.0) / (n_permut + 1.0)
            logger.info(f'Scoary1: done! {gene} {emp_p}')

        pvals.append(emp_p)
    return pvals
