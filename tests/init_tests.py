import json
import numpy as np
import pandas as pd
from os.path import dirname, exists
from scipy.spatial import distance
from scipy.stats import fisher_exact, boschloo_exact

from unittest import TestCase

# set up logging
import logging

logging.basicConfig()
# logging.getLogger().setLevel(logging.INFO)

ROOT = dirname(dirname(__file__))

DATA = {
    'generated': {
        'genes': 'Gene_presence_absence.csv',
        'traits': 'Trait.csv',
    },
    'tetracycline': {
        'genes': 'Gene_presence_absence.csv',
        'gene-info': 'gene-info.tsv',
        'traits': 'Tetracycline_resistance.csv',
        'traits-numeric': 'Tetracycline_resistance_numeric.csv',
        'restrict_to': 'Restrict_to.csv',
        'tree': 'ExampleTree.nwk',
        'treelist': 'expected_result.json',
        'tdm': 'tetracycline_TDM.csv',
        'scoary1-result': 'best_fisher_permute100.results.csv',
    },
    'small_ds': {
        'genes': 'pres_abs.csv',
        'traits': 'trait.csv',
        't1': 't1.results.csv',
        't2': 't2.results.csv',
    },
    'bigger_ds': {
        'traits': 'trait_trees.csv',
        'genes': 'pres_abs.csv',
        'tree': 'newick.nwk',
        'result-t1': 't1.results.csv',
        'result-t2': 't2.results.csv',
    },
    'new_ds': {
        'genes-og': 'Orthogroups.tsv',
        'genes-hog': 'N0.tsv',
        'genes-hog-info': 'N0_best_names.tsv',
        'isolate-meta': 'isolate-meta.tsv',
        'traits-lc-binary': 'LC-binary.tsv',
        'traits-lc': 'LC.tsv',
        'traits-lc-meta': 'LC-meta.tsv',
        'traits-gc-vol': 'GC-VOL.tsv',
        'traits-gc-vol-meta': 'none-meta.tsv',
    },
    'full_ds': {
        'genes': 'N0.tsv',
        'gene-info': 'N0_best_names.tsv',
        'isolate-info': 'isolate_info.tsv',
        'traits': 'traits.tsv',
        'trait-info': 'trait_info.tsv',
    },

}

VCF_DATA = {
    'vcf': {
        'vcf': 'Example.vcf',
        'traits': 'ExampleVCFTrait.csv'
    }
}

roary_ignore = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate',
                'Genome fragment', 'Order within fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
orthofinder_ignore = ['OG', 'Gene Tree Parent Clade']


def get_path(ds: str, key: str, data=DATA):
    return f'{ROOT}/data/{ds}/{data[ds][key]}'


def get_json(ds: str, key: str, data=DATA):
    with open(get_path(ds, key)) as f:
        return json.load(f)


def is_equivalent(a, b):
    if np.isinf(a) and np.isinf(b):
        return True
    if np.isnan(a) and np.isnan(b):
        return True
    return np.isclose(a, b)


def is_equivalent_tree(a, b) -> bool:
    if type(a) is str or type(b) is str:
        return a == b
    else:
        return (
                       is_equivalent_tree(a[0], b[0]) and is_equivalent_tree(a[1], b[1])
               ) or (
                       is_equivalent_tree(a[0], b[1]) and is_equivalent_tree(a[1], b[0])
               )


def get_tempdir_path() -> str:
    # template = '/tmp/scoary-test-outdir-{i}'
    # i = 0
    # while exists(template.format(i=i)):
    #     i += 1
    #
    # tempdir_path = template.format(i=i)

    tempdir_path = '/home/thomas/PycharmProjects/scoary-2/TEST_OUTPUT'

    logging.warning(f'Using this tempdir: file://{tempdir_path}')

    return tempdir_path
