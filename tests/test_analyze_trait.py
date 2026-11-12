import pandas as pd
import pytest
from scipy.stats import fisher_exact

from scoary.analyze_trait import init_result_df, create_test_df, add_odds_ratio, multiple_testing_correction, \
    pair_picking
from scoary.ScoaryTree import ScoaryTree
from conftest import TEST_DATA


@pytest.fixture
def small():
    #            s1     s2     s3     s4     s5     s6
    genes = pd.DataFrame({
        'pos': [True, True, True, False, False, False],  # perfect positive association
        'neg': [False, False, False, True, True, True],  # perfect negative association
        'mix': [True, False, True, True, False, False],
        'core': [True] * 6,
    }, index=['s1', 's2', 's3', 's4', 's5', 's6']).T
    genes.index.name = 'Gene'  # as set by load_genes
    trait = pd.Series([True, True, True, False, False, False],
                      index=['s1', 's2', 's3', 's4', 's5', 's6'], dtype='boolean')
    return genes, trait


def test_init_result_df_contingency(small):
    genes, trait = small
    df = init_result_df(genes, trait).set_index('Gene')

    assert 'core' not in df.index, 'genes present in all isolates must be dropped'
    assert df.loc['pos', '__contingency_table__'] == (3, 0, 0, 3)
    assert df.loc['neg', '__contingency_table__'] == (0, 3, 3, 0)
    assert df.loc['mix', '__contingency_table__'] == (2, 1, 1, 2)
    for col in ['g+t+', 'g+t-', 'g-t+', 'g-t-']:
        assert (df[col] >= 0).all()
    assert (df[['g+t+', 'g+t-', 'g-t+', 'g-t-']].sum(axis=1) == 6).all()


def test_sensitivity_specificity(small):
    genes, trait = small
    df = init_result_df(genes, trait).set_index('Gene')
    assert df.loc['pos', ['sensitivity', 'specificity']].tolist() == [100, 100]
    # regression test for #7: negatively correlated genes are scored by their absence
    assert df.loc['neg', ['sensitivity', 'specificity']].tolist() == [100, 100]
    assert df.loc['mix', 'sensitivity'] == pytest.approx(200 / 3)


def test_init_result_df_requires_boolean(small):
    genes, trait = small
    with pytest.raises(AssertionError):
        init_result_df(genes, trait.astype(bool))


def test_create_test_df(small):
    genes, trait = small
    result_df = init_result_df(genes, trait)
    test_df = create_test_df(result_df)
    assert test_df['fisher_p'].is_monotonic_increasing
    for _, row in test_df.iterrows():
        a, b, c, d = row['__contingency_table__']
        assert row['fisher_p'] == pytest.approx(fisher_exact([[a, b], [c, d]])[1])


def test_add_odds_ratio():
    test_df = pd.DataFrame({'__contingency_table__': [(4, 1, 2, 8)]})
    add_odds_ratio(test_df)
    assert test_df['odds_ratio'].iloc[0] == pytest.approx((4 * 8) / (1 * 2))


def test_multiple_testing_correction_native():
    df = pd.DataFrame({'p': [0.001, 0.01, 0.04, 0.2]})
    out = multiple_testing_correction(df.copy(), 'p', 'q', method='native', cutoff=0.05)
    assert out['p'].tolist() == [0.001, 0.01, 0.04]  # native: filter on raw p-value
    assert out['q'].tolist() == pytest.approx([0.004, 0.04, 0.16])  # q = bonferroni


def test_multiple_testing_correction_bonferroni():
    df = pd.DataFrame({'p': [0.001, 0.01, 0.04, 0.2]})
    out = multiple_testing_correction(df.copy(), 'p', 'q', method='bonferroni', cutoff=0.05)
    assert out['p'].tolist() == [0.001, 0.01]


def test_tetracycline_matches_scoary_1(tetracycline_genes, tetracycline_traits, tetracycline_tree_list):
    """Regression test: contingency tables, odds ratios and pair picking must match Scoary 1."""
    trait = tetracycline_traits['Tetracycline_resistance']
    result_df = init_result_df(tetracycline_genes, trait_series=trait.astype('boolean'))
    result_df = add_odds_ratio(result_df)
    tree = ScoaryTree.from_list(tetracycline_tree_list)
    result_df = pair_picking(result_df, tetracycline_genes, tree=tree, label_to_trait=trait)
    result_df.set_index('Gene', inplace=True)

    expected = pd.read_csv(TEST_DATA / 'fisher_permute100.results.csv')
    assert len(expected) > 0
    for _, row in expected.iterrows():
        new = result_df.loc[row.Gene]
        assert (row.Number_pos_present_in, row.Number_neg_present_in,
                row.Number_pos_not_present_in, row.Number_neg_not_present_in) == \
               tuple(int(new[c]) for c in ('g+t+', 'g+t-', 'g-t+', 'g-t-')), row.Gene
        assert new.odds_ratio == pytest.approx(row.Odds_ratio), row.Gene
        # Scoary 1 always reports sensitivity/specificity for presence of the gene. Scoary 2 (since #7)
        # reports them for absence if the gene is negatively associated, i.e. 100 - Scoary 1's value.
        if row.Sensitivity + row.Specificity >= 100:
            exp_sens, exp_spec = row.Sensitivity, row.Specificity
        else:
            exp_sens, exp_spec = 100 - row.Sensitivity, 100 - row.Specificity
        assert new.sensitivity == pytest.approx(exp_sens), row.Gene
        assert new.specificity == pytest.approx(exp_spec), row.Gene
        assert new.contrasting == row.Max_Pairwise_comparisons, row.Gene
        assert new.supporting == row.Max_supporting_pairs, row.Gene
        assert new.opposing == row.Max_opposing_pairs, row.Gene
        assert new.best == pytest.approx(row.Best_pairwise_comp_p), row.Gene
        assert new.worst == pytest.approx(row.Worst_pairwise_comp_p), row.Gene
