"""End-to-end tests of the full pipeline on the (tracked) tetracycline example."""
import pandas as pd
import pytest

from scoary.scoary import scoary, find_duplicates
from conftest import TETRACYCLINE, TEST_DATA


def test_find_duplicates():
    traits = pd.DataFrame({
        'a': [True, False, True],
        'b': [False, False, True],
        'a_copy': [True, False, True],
        'a_copy2': [True, False, True],
        'c': [False, True, True],
    })
    assert find_duplicates(traits).to_dict() == {'a_copy': 'a', 'a_copy2': 'a'}


def test_find_duplicates_none():
    traits = pd.DataFrame({'a': [True, False], 'b': [False, True]})
    assert len(find_duplicates(traits)) == 0


@pytest.mark.slow
@pytest.mark.parametrize('n_cpus', [1, 2])
def test_tetracycline_end_to_end(tmp_path, n_cpus):
    outdir = tmp_path / 'out'
    scoary(
        genes=str(TETRACYCLINE / 'Gene_presence_absence.csv'),
        traits=str(TETRACYCLINE / 'Tetracycline_resistance.csv'),
        gene_info=str(TEST_DATA / 'gene-info.tsv'),
        multiple_testing='bonferroni:0.05',
        n_permut=200,
        random_state=42,
        n_cpus=n_cpus,
        outdir=str(outdir),
    )
    summary = pd.read_csv(outdir / 'summary.tsv', sep='\t', index_col=0)
    assert 'Tetracycline_resistance' in summary.index

    result = pd.read_csv(outdir / 'traits' / 'Tetracycline_resistance' / 'result.tsv', sep='\t')
    assert len(result) > 0
    assert (result['fisher_q'] <= 0.05).all()
    assert result['fisher_p'].is_monotonic_increasing

    # the known resistance-associated genes of the example dataset should top the list
    assert result['Gene'].iloc[0] in {'TetRCG', 'tetM', 'tetO', 'tetR'} or result['odds_ratio'].iloc[0] > 10
