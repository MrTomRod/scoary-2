# Scoary 2

Associate genes with traits!

# Installation

1) [Python](https://pypi.org/project/scoary-2/) 3.10+: `pip install scoary-2`
2) [Docker](https://hub.docker.com/r/troder/scoary-2): `docker pull troder/scoary-2`

# Usage

Examples:

1) Dataset from Scoary 1

```shell
# Dataset from Scoary 1: genes in Roary gene count format
scoary \
    --genes Gene_presence_absence.csv \
    --traits Tetracycline_resistance.csv \
    --outdir out \
    --n-permut 1000
```

See output [here](https://gfv-oberburg.ch/FILES/SCOARY2_TETR/overview.html).

2) OrthoFinder-based data

```shell
 scoary \
    --genes N0.tsv \
    --gene-data-type 'gene-list:\t' \
    --gene-info N0_best_names.tsv \
    --traits traits.tsv \
    --trait-data-type 'gaussian:\t' \
    --trait-info trait_info.tsv \
    --isolate-info isolate_info.tsv \
    --n-permut 200 \
    --n-cpus 7 \
    --outdir out
# the following are optional: gene-info, trait-info, isolate-info
```

See output of a huge, unfiltered dataset [here](https://gfv-oberburg.ch/FILES/SCOARY2.2/overview.html).
(Takes a few seconds to load.)

## Docker cookbook

Run interactive shell:

```shell
docker run  \
    -u $(id -u ${USER}):$(id -g ${USER})  \
    -it --rm  \
    -v /path/to/data:/data:Z  \
    troder/scoary-2  \
    /bin/bash
```

## Options

Below is the output of `scoary --help`:

```text
POSITIONAL ARGUMENTS
    GENES
        Type: str
        Path to gene presence/absence table: columns=isolates, rows=genes
    TRAITS
        Type: str
        Path to trait presence/absence table: columns=traits, rows=isolates
    OUTDIR
        Type: str
        Directory to place output files

FLAGS
    --multiple_testing_fisher=MULTIPLE_TESTING_FISHER
        Type: str
        Default: 'bonferroni:0.999'
        "method:cutoff" for filtering genes after Fisher's test, where cutoff is a number and method is one of [bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by,  fdr_tsbh, fdr_tsbky]
    --multiple_testing_picking=MULTIPLE_TESTING_PICKING
        Type: str
        Default: 'bonferroni:0.999'
        "method:cutoff" for filtering genes after the pairwise comparisons algorithm, where cutoff is a number and method is one of [bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by,  fdr_tsbh, fdr_tsbky]
    --gene_info=GENE_INFO
        Type: Optional[str]
        Default: None
        Path to file that describes genes: columns=arbitrary properties, rows=genes
    --trait_info=TRAIT_INFO
        Type: Optional[str]
        Default: None
        Path to file that describes traits: columns=arbitrary properties, rows=traits
    --isolate_info=ISOLATE_INFO
        Type: Optional[str]
        Default: None
        Path to file that describes isolates: columns=arbitrary properties, rows=isolates
    --newicktree=NEWICKTREE
        Type: Optional[str]
        Default: None
        Path to a custom tree in Newick format
    --pairwise=PAIRWISE
        Type: bool
        Default: True
        If False, only perform Fisher's test. If True, also perform pairwise comparisons algorithm.
    --n_permut=N_PERMUT
        Type: int
        Default: 500
        Post-hoc label-switching test: perform N permutations of the phenotype by random label switching. Low p-values suggest that the effect is not merely lineage-specific.
    --restrict_to=RESTRICT_TO
        Type: Optional[str]
        Default: None
        Comma-separated list of isolates to which to restrict this analysis
    --ignore=IGNORE
        Type: Optional[str]
        Default: None
        Comma-separated list of isolates to be ignored for this analysis
    --n_cpus=N_CPUS
        Type: int
        Default: 1
    --trait_data_type=TRAIT_DATA_TYPE
        Type: str
        Default: 'binary:,'
        "<method>:<?cutoff>:<?covariance_type>:<?alternative>:<?delimiter>" How to read the traits table. Example: "gene-list:\t" for OrthoFinder N0.tsv table
    --gene_data_type=GENE_DATA_TYPE
        Type: str
        Default: 'gene-count:,'
        "<data_type>:<?delimiter>" How to read the genes table. Example: "gene-list:\t" for OrthoFinder N0.tsv table
    --random_state=RANDOM_STATE
        Type: Optional[int]
        Default: None
        Set a fixed seed for the random number generator
    --limit_traits=LIMIT_TRAITS
        Type: Optional[int, int]
        Default: None
        Limit the analysis to traits n to m. Useful for debugging. Example: "(0, 10)"
```

# Overview of the algorithm

![algorithm flowchart](media/ScoaryWorkflow.drawio.svg)

# Usage of the picking algorithm in Python

Python bindings to the _pairwise comparisons algorithm_, as described in 
[Read, 1995](https://doi.org/10.1006/jtbi.1995.0047), 
[Maddison, 2000](https://doi.org/10.1006/jtbi.1999.1050) and 
[Brynildsrud, 2016](https://doi.org/10.1186/s13059-016-1108-8).

<details>

  <summary>Simple pair picking</summary>

```python
from pprint import pprint
from scoary import ScoaryTree, pick_single, print_tree

tree = [['isolate1', 'isolate2'], [['isolate3', 'isolate4'], ['isolate5', 'isolate6']]]

label_to_trait_a = {
    'isolate1': True,
    'isolate2': False,
    'isolate3': True,
    'isolate4': False,
    'isolate5': True,
    'isolate6': False,
}

label_to_trait_b = {
    'isolate1': True,
    'isolate2': False,
    'isolate3': True,
    'isolate4': False,
    'isolate5': True,
    'isolate6': False,
}

print_tree(
    ScoaryTree.from_list(tree),
    label_to_trait_a, label_to_trait_b
)
#       /-11_isolate1
#    /-|
#   |   \-00_isolate2
#   |
# --|      /-11_isolate3
#   |   /-|
#   |  |   \-00_isolate4
#    \-|
#      |   /-11_isolate5
#       \-|
#          \-00_isolate6

result = pick_single(tree, label_to_trait_a, label_to_trait_b, calc_pvals=True)
pprint(result)
# {'best_pval': 0.25,
#  'max_contrasting_pairs': 3,
#  'max_opposing_pairs': 0,
#  'max_supporting_pairs': 3,
#  'worst_pval': 0.25}
 ```

</details>

<details>

  <summary>Parallel pair picking</summary>

This takes advantage of [Numba](https://numba.pydata.org/) optimizations.

```python
import pandas as pd
from scoary import pick

tree = [['isolate1', 'isolate2'], ['isolate3', 'isolate4']]

# e.g. phenotype
label_to_trait_a = {
    'isolate1': True,
    'isolate2': False,
    'isolate3': False,
    'isolate4': True,
}

# e.g. presence/absence of genes
trait_b_df = pd.DataFrame(
    columns=['isolate1', 'isolate2', 'isolate3', 'isolate4'],
    data=[
        [True, True, False, False],  # gene 1
        [True, False, True, False],  # gene 2
        [True, False, False, True],  # ...
        [False, True, True, False],
        [False, True, False, True],
        [False, True, False, True],
        [False, True, False, True],
        [False, True, False, True],
    ]
)

max_contr, max_suppo, max_oppos, best, worst = pick(
    tree=tree,
    label_to_trait_a=label_to_trait_a,
    trait_b_df=trait_b_df,
    calc_pvals=True
)

print(f'{max_contr=}\n{max_suppo=}\n{max_oppos=}\n{best=}\n{worst=}')
# max_contr=array([1, 2, 2, 2, 2, 2, 2, 2])
# max_suppo=array([1, 1, 2, 0, 1, 1, 1, 1])
# max_oppos=array([1, 1, 0, 2, 1, 1, 1, 1])
# best=array([1. , 1. , 0.5, 0.5, 1. , 1. , 1. , 1. ])
# worst=array([1. , 1. , 0.5, 0.5, 1. , 1. , 1. , 1. ])
```

</details>

# Todo:

- [X] Binarize traits during multiprocessing
- [ ] Idea: min_qval * min_pval_empirical as new score
- [ ] Improve messaging and logging
- [ ] Add links to `summary.tsv` / `coverage-matrix.tsv` / `result.tsv`
- [ ] Benchmark?
- [ ] Refactoring
- [X] `overview.html`: go to function, popover always left
- [ ] `overview.html` & `trait.html`: set proper title
- [ ] `overview.html` & `trait.html`: link to table download
- [ ] Improve `README.md`
