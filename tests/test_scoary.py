from init_tests import *

from scoary.scoary import *


class TestScoary(TestCase):
    def test_scoary(self):
        scoary(
            genes=get_path('tetracycline', 'genes'),
            traits=get_path('tetracycline', 'traits'),
            start_col=13,
        )

    def test_load_genes(self):
        res = load_genes(get_path('generated', 'genes'), delimiter=',', start_col=0)
        print(res)
        res = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        print(res)

    def test_load_traits(self):
        res = load_traits(get_path('generated', 'traits'), delimiter=',')
        print(res)
        res = load_traits(get_path('tetracycline', 'traits'), delimiter=',')
        print(res)

    def test_same_hemming_result(self):
        """
        Check if old scoary generates the same data (hamming similarity matrix)
        """
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        tdm_new = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')))
        tdm_old = np.flip(pd.read_csv(get_path('tetracycline', 'tdm'), index_col=0).values)  # has to be flipped
        np.fill_diagonal(tdm_old, 0)  # diagonal should be 0, not 1
        tdm_old = pd.DataFrame(tdm_old)
        self.assertTrue(np.isclose(tdm_old, tdm_new).all())

    def test_create_result_df(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        result_df = init_result_df(genes_df, trait_pos=set(genes_df.columns[:10]), trait_neg=set(genes_df.columns[90:]))
        self.assertEqual(result_df.columns.tolist(), ['c1r1', 'c2r1', 'c1r2', 'c2r2', '__contingency_table__'])

    def test_contingency_test(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        result_df = init_result_df(genes_df, trait_pos=set(genes_df.columns[:11]), trait_neg=set(genes_df.columns[89:]))
        test_df = create_test_df(result_df=result_df)
        self.assertEqual(['__contingency_table__', 'pval'], test_df.columns.tolist())
        print(f"Done: minpval={test_df.pval.min()}")

    def test_odds_ratio(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        genes_df = genes_df[:100]  # only first 100 rows
        test_df = init_result_df(genes_df, trait_pos=set(genes_df.columns[:11]), trait_neg=set(genes_df.columns[89:]))

        # apply function
        test_df = add_odds_ratio(test_df)
        self.assertEqual(test_df.columns.tolist(), ['c1r1', 'c2r1', 'c1r2', 'c2r2', '__contingency_table__', 'odds_ratio'])

        # calculate odds_ratio with fisher_exact
        fisher_ors = test_df.apply(lambda row: fisher_exact([[row['c1r1'], row['c2r1']], [row['c1r2'], row['c2r2']]])[0], axis=1)

        # check if result is identical
        for manual_or, fisher_or in zip(test_df['odds_ratio'], fisher_ors):
            self.assertTrue(is_equivalent(manual_or, fisher_or))
        print(test_df.columns)

    def test_odds_ratio_2(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        traits_df = load_traits(get_path('tetracycline', 'traits'), delimiter=',')
        trait_df = traits_df['Tetracycline_resistance']
        test_df = init_result_df(genes_df, trait_pos=set(trait_df[trait_df == 1].index), trait_neg=set(trait_df[trait_df == 0].index))

        # apply function
        test_df = add_odds_ratio(test_df)
        test_data = test_df.odds_ratio.to_dict()

        # load expected result from scoary 1
        expected_result = pd.read_csv(get_path('tetracycline', 'scoary1-result'))

        # check if result is identical
        for i, row in expected_result.iterrows():
            table = (row.Number_pos_present_in,
                     row.Number_neg_present_in,
                     row.Number_pos_not_present_in,
                     row.Number_neg_not_present_in)
            new_row = test_df.loc[row.Gene]
            new_table = tuple(int(new_row[c]) for c in ('c1r1', 'c2r1', 'c1r2', 'c2r2'))

            self.assertEqual(table, new_table)
            self.assertTrue(np.isclose(row.Odds_ratio, test_data[row.Gene]),
                            msg=f'Failed to calculate odds_ratio for {row.Gene}: {row.Odds_ratio} != {test_data[row.Gene]}')
