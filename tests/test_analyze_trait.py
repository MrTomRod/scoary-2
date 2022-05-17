from init_tests import *
from scoary.utils import get_label_to_trait, get_all_label_to_gene
from scoary.ScoaryTree import ScoaryTree
from scoary.load_genes import load_genes
from scoary.load_traits import load_traits
from scoary.analyze_trait import init_result_df, create_test_df, add_odds_ratio, pair_picking


def generate_fake_traits(genes_df: pd.DataFrame) -> {str: bool}:
    label_to_trait = {}
    label_to_trait.update({l: True for l in genes_df.columns[:11]})
    label_to_trait.update({l: False for l in genes_df.columns[89:]})
    return label_to_trait


class TestScoary(TestCase):
    def test_create_result_df(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=roary_ignore)
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        self.assertEqual(
            result_df.columns.tolist(),
            ['Gene', 'g+t+', 'g+t-', 'g-t+', 'g-t-', '__contingency_table__', 'sensitivity', 'specificity']
        )

    def test_contingency_test(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=roary_ignore)
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        test_df = create_test_df(result_df=result_df)
        self.assertEqual(['__contingency_table__', 'pval'], test_df.columns.tolist())
        print(f"Done: minpval={test_df.pval.min()}")

    def test_odds_ratio(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=roary_ignore)
        genes_df = genes_df[:100]  # only first 100 rows
        test_df = init_result_df(genes_df, generate_fake_traits(genes_df))

        # apply function
        test_df = add_odds_ratio(test_df)
        self.assertEqual(
            test_df.columns.tolist(),
            ['Gene', 'g+t+', 'g+t-', 'g-t+', 'g-t-', '__contingency_table__', 'sensitivity', 'specificity',
             'odds_ratio']
        )

        # calculate odds_ratio with fisher_exact
        fisher_ors = test_df.apply(
            lambda row: fisher_exact([[row['g+t+'], row['g+t-']], [row['g-t+'], row['g-t-']]])[0], axis=1)

        # check if result is identical
        for manual_or, fisher_or in zip(test_df['odds_ratio'], fisher_ors):
            self.assertTrue(is_equivalent(manual_or, fisher_or))

    def test_tetracycline(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=roary_ignore)
        _, traits_df = load_traits(get_path('tetracycline', 'traits'), trait_data_type='binary:,')
        trait_series = traits_df['Tetracycline_resistance']

        # calculate sensitivity and specificity
        test_df = init_result_df(genes_df, label_to_trait={l: bool(v) for l, v in trait_series.items() if v in (0, 1)})
        # calculate odds_ratio
        test_df = add_odds_ratio(test_df)
        # calculate pairwise comparisons
        all_label_to_gene = get_all_label_to_gene(genes_df)
        tree = ScoaryTree.from_list(get_json('tetracycline', 'treelist')['as_list'])
        assert set(tree.labels()) == set(all_label_to_gene['TetRCG'])
        label_to_trait = get_label_to_trait(trait_series)
        test_df = pair_picking(test_df, genes_df, tree=tree, label_to_trait=label_to_trait)

        # load expected result from scoary 1
        expected_result = pd.read_csv(get_path('tetracycline', 'scoary1-result'))

        test_df.set_index('Gene', inplace=True)

        # check if result is identical
        for i, row in expected_result.iterrows():
            table = (row.Number_pos_present_in,
                     row.Number_neg_present_in,
                     row.Number_pos_not_present_in,
                     row.Number_neg_not_present_in)
            new_row = test_df.loc[row.Gene]
            new_table = tuple(int(new_row[c]) for c in ('g+t+', 'g+t-', 'g-t+', 'g-t-'))

            self.assertEqual(table, new_table)
            self.assertAlmostEqual(
                row.Odds_ratio, new_row.odds_ratio,
                msg=f'Failed to calculate odds_ratio for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}'
            )
            self.assertAlmostEqual(
                row.Sensitivity, new_row.sensitivity,
                msg=f'Failed to calculate sensitivity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}'
            )
            self.assertAlmostEqual(
                row.Specificity, new_row.specificity,
                msg=f'Failed to calculate specificity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}'
            )

            xx = [
                (row.Max_Pairwise_comparisons, new_row.contrasting),
                (row.Max_supporting_pairs, new_row.supporting),
                (row.Max_opposing_pairs, new_row.opposing),
                (row.Best_pairwise_comp_p, new_row.best),
                (row.Worst_pairwise_comp_p, new_row.worst)
            ]
            try:
                self.assertEqual(row.Max_Pairwise_comparisons, new_row.contrasting)
                self.assertEqual(row.Max_supporting_pairs, new_row.supporting)
                self.assertEqual(row.Max_opposing_pairs, new_row.opposing)
                self.assertAlmostEqual(row.Best_pairwise_comp_p, new_row.best)
                self.assertAlmostEqual(row.Worst_pairwise_comp_p, new_row.worst)
            except Exception as e:
                print(i, row.Gene, xx)
                self.fail(msg=str(e))
