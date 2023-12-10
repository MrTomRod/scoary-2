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

roary_ignore = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate',
                'Genome fragment', 'Order within fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
orthofinder_ignore = ['OG', 'Gene Tree Parent Clade']


def get_json(path: str):
    with open(path) as f:
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
