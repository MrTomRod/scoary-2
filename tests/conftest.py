"""
Shared fixtures for the scoary-2 unit tests.

Only data that is tracked in git is used here (``data/tetracycline`` and ``tests/data``),
so the suite runs on a fresh clone, from any working directory.
"""
import json
import os
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TETRACYCLINE = REPO / 'data' / 'tetracycline'
TEST_DATA = Path(__file__).resolve().parent / 'data'

# columns of the Roary output that are metadata, not isolates
ROARY_IGNORE = [
    'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate',
    'Genome fragment', 'Order within fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
    'Min group size nuc', 'Max group size nuc', 'Avg group size nuc'
]

# don't copy the mgwas-data-exploration-app into every test output
os.environ['MGWAS_LINK_ONLY'] = 'true'


@pytest.fixture(scope='session')
def tetracycline_tree_list() -> list:
    with open(TETRACYCLINE / 'expected_result.json') as f:
        return json.load(f)['as_list']


@pytest.fixture(scope='session')
def tetracycline_genes():
    from scoary.load_genes import load_genes
    _, genes_df = load_genes(
        str(TETRACYCLINE / 'Gene_presence_absence.csv'), gene_data_type='gene-count', ignore=ROARY_IGNORE
    )
    return genes_df


@pytest.fixture(scope='session')
def tetracycline_traits():
    from scoary.load_traits import load_traits
    _, traits_df = load_traits(str(TETRACYCLINE / 'Tetracycline_resistance.csv'), trait_data_type='binary:,')
    return traits_df


def is_equivalent_tree(a, b) -> bool:
    """Compare two list-trees, ignoring the order of children."""
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    return (is_equivalent_tree(a[0], b[0]) and is_equivalent_tree(a[1], b[1])) or \
        (is_equivalent_tree(a[0], b[1]) and is_equivalent_tree(a[1], b[0]))
