import itertools
import logging

import pytest
from fast_fisher.fast_fisher_cython import test1t as fisher_exact_two_tailed
from scipy.stats import fisher_exact

from scoary.utils import parse_correction, split_into_parts, fisher_id, is_int, is_float, load_info_file
from conftest import TEST_DATA


@pytest.mark.parametrize('s, expected', [
    ('native:0.05', ('native', 0.05)),
    ('fdr_bh:0.3', ('fdr_bh', 0.3)),
    ('bonferroni', ('bonferroni', float('inf'))),
])
def test_parse_correction(s, expected):
    assert parse_correction(s, 'mt') == expected


@pytest.mark.parametrize('bad', ['not-a-method:0.05', 'native:abc'])
def test_parse_correction_rejects(bad):
    with pytest.raises(AssertionError):
        parse_correction(bad, 'mt')


@pytest.mark.parametrize('n_items, n_parts', [(10, 3), (3, 5), (0, 2), (7, 7), (100, 1)])
def test_split_into_parts(n_items, n_parts):
    items = list(range(n_items))
    parts = split_into_parts(items, n_parts)
    assert len(parts) == n_parts
    assert list(itertools.chain(*parts)) == items  # order preserved, nothing lost
    sizes = [len(p) for p in parts]
    assert max(sizes) - min(sizes) <= 1  # balanced


def test_is_int_is_float():
    assert is_int('3') and not is_int('3.5') and not is_int('x')
    assert is_float('3') and is_float('3.5') and is_float('1e-5') and not is_float('x')


@pytest.mark.parametrize('table', [(5, 1, 2, 8), (0, 10, 10, 0), (3, 3, 3, 3), (12, 0, 1, 7)])
def test_fisher_id_is_invariant(table):
    """All eight symmetric variants of a 2x2 table share one id and one p-value."""
    a, b, c, d = table
    variants = [(a, b, c, d), (a, c, b, d), (b, a, d, c), (b, d, a, c),
                (c, a, d, b), (c, d, a, b), (d, b, c, a), (d, c, b, a)]
    ids = {fisher_id(*v) for v in variants}
    assert len(ids) == 1
    pvals = [fisher_exact([[v[0], v[1]], [v[2], v[3]]])[1] for v in variants]
    assert pvals == pytest.approx([pvals[0]] * 8)


@pytest.mark.parametrize('table', [(5, 1, 2, 8), (0, 10, 10, 0), (3, 3, 3, 3), (20, 4, 7, 30)])
def test_fast_fisher_matches_scipy(table):
    a, b, c, d = table
    assert fisher_exact_two_tailed(a, b, c, d) == pytest.approx(fisher_exact([[a, b], [c, d]])[1])


def test_load_info_file():
    df = load_info_file(logging.getLogger(), str(TEST_DATA / 'gene-info.tsv'), merge_col='Gene')
    assert df.index.name == 'Gene'
    assert df.index.is_unique


def test_load_info_file_wrong_header():
    with pytest.raises(AssertionError, match='first column must be named'):
        load_info_file(logging.getLogger(), str(TEST_DATA / 'gene-info.tsv'), merge_col='Isolate')
