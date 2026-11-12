import pytest

from scoary.load_genes import load_genes, parse_params
from conftest import TETRACYCLINE, ROARY_IGNORE

COUNT_CSV = """\
Gene,s1,s2,s3,s4
core,1,1,1,1
absent,0,0,0,0
g1,1,0,2,0
g2,0,3,0,1
"""

LIST_TSV = """\
OG\ts1\ts2\ts3\ts4
core\ta1\tb1\tc1\td1
g1\ta2\t\tc2, c3\t
g2\t\tb2\t\td2
"""


@pytest.fixture
def count_file(tmp_path):
    p = tmp_path / 'genes.csv'
    p.write_text(COUNT_CSV)
    return str(p)


@pytest.fixture
def list_file(tmp_path):
    p = tmp_path / 'genes.tsv'
    p.write_text(LIST_TSV)
    return str(p)


@pytest.mark.parametrize('s, expected', [
    ('gene-count', ('gene-count', ',')),
    ('gene-count:;', ('gene-count', ';')),
    ('gene-list:\t', ('gene-list', '\t')),
    ('GENE-LIST:,', ('gene-list', ',')),
])
def test_parse_params(s, expected):
    assert parse_params(s) == expected


@pytest.mark.parametrize('bad', ['gene-foo', 'gene-count:,:x'])
def test_parse_params_rejects(bad):
    with pytest.raises(AssertionError):
        parse_params(bad)


def test_load_gene_count(count_file):
    orig_df, bool_df = load_genes(count_file, 'gene-count')
    # core and absent genes are dropped, counts >= 1 become True
    assert bool_df.index.tolist() == ['g1', 'g2']
    assert bool_df.loc['g1'].tolist() == [True, False, True, False]
    assert bool_df.loc['g2'].tolist() == [False, True, False, True]
    assert bool_df.index.name == orig_df.index.name == 'Gene'
    # original counts are kept
    assert orig_df.loc['g2', 's2'] == 3


def test_load_gene_list(list_file):
    _, bool_df = load_genes(list_file, 'gene-list:\t')
    assert bool_df.index.tolist() == ['g1', 'g2']
    assert bool_df.loc['g1'].tolist() == [True, False, True, False]
    assert bool_df.loc['g2'].tolist() == [False, True, False, True]


def test_restrict_to_and_ignore(count_file):
    _, bool_df = load_genes(count_file, 'gene-count', restrict_to=['s1', 's2', 's3'])
    assert bool_df.columns.tolist() == ['s1', 's2', 's3']

    _, bool_df = load_genes(count_file, 'gene-count', ignore=['s4'])
    assert bool_df.columns.tolist() == ['s1', 's2', 's3']


@pytest.mark.parametrize('kwargs', [{'restrict_to': ['s1', 'nope']}, {'ignore': ['nope']}])
def test_unknown_strains_fail(count_file, kwargs):
    with pytest.raises(AssertionError):
        load_genes(count_file, 'gene-count', **kwargs)


def test_nan_in_count_file_fails(tmp_path):
    p = tmp_path / 'genes.csv'
    p.write_text('Gene,s1,s2\ng1,1,\n')
    with pytest.raises(AssertionError, match='NaN'):
        load_genes(str(p), 'gene-count')


def test_duplicate_genes_fail(tmp_path):
    p = tmp_path / 'genes.csv'
    p.write_text('Gene,s1,s2\ng1,1,0\ng1,0,1\n')
    with pytest.raises(AssertionError, match='index not unique'):
        load_genes(str(p), 'gene-count')


def test_tetracycline(tetracycline_genes):
    assert tetracycline_genes.shape[1] == 100
    assert tetracycline_genes.dtypes.eq(bool).all()
    row_sums = tetracycline_genes.sum(axis=1)
    assert ((row_sums > 0) & (row_sums < 100)).all()
