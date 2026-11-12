import numpy as np
import pandas as pd
import pytest

from scoary.load_traits import load_traits, parse_params
from conftest import TEST_DATA

BINARY_CSV = """\
,t1,t2
s1,1,0
s2,0,1
s3,1,
s4,0,1
"""


@pytest.mark.parametrize('s, expected', [
    ('binary', ('binary', ',', 0.85, 'tied', 'skip')),
    ('binary:\t', ('binary', '\t', 0.85, 'tied', 'skip')),
    ('gaussian:kmeans:.9:full:;', ('gaussian', ';', 0.9, 'full', 'kmeans')),
    ('kmeans', ('kmeans', ',', 0.85, 'tied', 'skip')),
])
def test_parse_params(s, expected):
    assert parse_params(s) == expected


@pytest.mark.parametrize('bad', ['foo', 'gaussian:0.3', 'gaussian:full:tied', 'binary:,:;'])
def test_parse_params_rejects(bad):
    with pytest.raises(AssertionError):
        parse_params(bad)


def test_load_binary(tmp_path):
    p = tmp_path / 'traits.csv'
    p.write_text(BINARY_CSV)
    numeric_df, traits_df = load_traits(str(p), 'binary:,')
    assert numeric_df is None
    assert traits_df.index.name == 'Isolate'
    assert traits_df.index.tolist() == ['s1', 's2', 's3', 's4']
    assert (traits_df.dtypes == 'boolean').all()
    assert traits_df['t1'].tolist() == [True, False, True, False]
    assert pd.isna(traits_df.loc['s3', 't2'])  # missing values are kept as NA


def test_load_binary_limit_and_filter(tmp_path):
    p = tmp_path / 'traits.csv'
    p.write_text(BINARY_CSV)
    _, traits_df = load_traits(str(p), 'binary:,', limit_traits=(1, 2))
    assert traits_df.columns.tolist() == ['t2']


def test_tetracycline(tetracycline_traits):
    assert tetracycline_traits.columns.tolist() == ['Tetracycline_resistance', 'Bogus_trait']
    assert len(tetracycline_traits) == 100


@pytest.mark.parametrize('method', ['kmeans', 'gaussian'])
def test_binarize_clearly_bimodal(tmp_path, method):
    rng = np.random.default_rng(0)
    low, high = rng.normal(0, .1, 20), rng.normal(10, .1, 20)
    df = pd.DataFrame({'trait': np.concatenate([low, high])}, index=[f's{i}' for i in range(40)])
    p = tmp_path / 'numeric.csv'
    df.to_csv(p)

    numeric_df, traits_df = load_traits(str(p), f'{method}:,', random_state=42)
    assert numeric_df.shape == (40, 1)
    binarized = traits_df['trait']
    # high values must be True, low values False
    assert binarized.iloc[20:].all()
    assert not binarized.iloc[:20].any()


def test_binarize_tetracycline_numeric():
    numeric_df, traits_df = load_traits(
        str(TEST_DATA / 'Tetracycline_resistance_numeric.csv'), 'gaussian:kmeans:,', random_state=42
    )
    assert traits_df.index.equals(numeric_df.index)
    for trait in traits_df.columns:
        t = traits_df[trait]
        values = numeric_df[trait]
        # every True isolate must have a larger value than every False isolate
        if t.notna().all() and t.any() and not t.all():
            assert values[t.astype(bool)].min() > values[~t.astype(bool)].max()
