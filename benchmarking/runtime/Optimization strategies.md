# Comparing different optimization strategies

I also tried original [Scoary's optimization](https://github.com/AdmiralenOla/Scoary/blob/b713e10fc1968488132f62652c6dba35636ca3e6/scoary/methods.py#L1360-L1363)
(breaking the permutations) instead of my approach, caching confidence intervals.

## Results

**Scoary:**

- break disabled: 41 minutes
- normal: 22 minutes

**Scoary2 (1 CPU):**

- cache disabled: 2:01
- normal: 1:12
- break instead of cache: 1:46

**Scoary2 (8 CPUs):**

- cache: 26 sec
- break: 39 sec

## Summary

My caching optimization appears to be better.

## Code

The code below replaces the permute_picking function in [permutations.py](/scoary/permutations.py)

Note: I have not thoroughly tested this code, so it may contain bugs.

```python
import scipy.stats as ss


def permute_picking(
        trait: str,
        tree: ScoaryTree,
        label_to_trait: pd.Series | dict,
        result_df: pd.DataFrame,
        genes_bool_df: pd.DataFrame,
        n_permut: int,
        random_state: int = None,
        batch_size: int = 50
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

        is_positively_correlated = row.supporting >= row.opposing
        estimator = (row.supporting if is_positively_correlated else row.opposing) / row.contrasting
        n_pos_assoc = n_pos if is_positively_correlated else n_neg

        r = 0
        for batch_start in range(0, n_permut, batch_size):
            batch_end = min(batch_start + batch_size, n_permut)
            batch_size_current = batch_end - batch_start
            # print(f"Processing {batch_start + 1}-{batch_end} ({batch_size_current} of {n_permut} items)")

            permuted_df = create_permuted_df(
                labels=labels, n_positive=n_pos_assoc,
                n_permut=batch_size_current, random_state=random_state
            )
            max_contr, max_suppo, max_oppos = pick(
                tree=tree.to_list, label_to_trait_a=label_to_gene,
                trait_b_df=permuted_df, calc_pvals=False
            )

            # Check how many estimators are higher than the unpermuted
            r += sum((max_suppo / max_contr) >= estimator)

            # If r indicates a p > 0.1 with a probability of 95%, abort
            if batch_end >= 30 and (1 - ss.binom.cdf(r, batch_end, 0.1)) < 0.05:
                pval = (r + 1) / (batch_end + 1)
                break

        else:
            pval = (r + 1) / (n_permut + 1)

        pvals.append(pval)

    return pvals
```