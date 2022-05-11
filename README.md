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
    --threads 7 \
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

# Overview of the algorithm

![algorithm flowchart](media/ScoaryWorkflow.drawio.svg)

# Todo:

- [X] Binarize traits during multiprocessing
- [ ] Idea: min_qval * min_pval_empirical as new score
- [ ] Improve messaging and logging
- [ ] Add links to `summary.tsv` / `coverage-matrix.tsv` / `result.tsv`
- [ ] Benchmark?
- [ ] Refactoring
- [ ] `overview.html` & `trait.html`: set proper title
- [ ] Improve `README.md`
