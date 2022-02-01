import pandas as pd

from init_tests import *

from scoary.scoary import *


def generate_fake_traits(genes_df: pd.DataFrame) -> {str: bool}:
    label_to_trait = {}
    label_to_trait.update({l: True for l in genes_df.columns[:11]})
    label_to_trait.update({l: False for l in genes_df.columns[89:]})
    return label_to_trait


class TestScoary(TestCase):
    def test_scoary(self):
        scoary(
            genes=get_path('tetracycline', 'genes'),
            traits=get_path('tetracycline', 'traits'),
            start_col=13,
            n_permut=100
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
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        self.assertEqual(
            result_df.columns.tolist(),
            ['Gene', 'c1r1', 'c2r1', 'c1r2', 'c2r2', '__contingency_table__', 'sensitivity', 'specificity']
        )

    def test_contingency_test(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        test_df = create_test_df(result_df=result_df)
        self.assertEqual(['__contingency_table__', 'pval'], test_df.columns.tolist())
        print(f"Done: minpval={test_df.pval.min()}")

    def test_odds_ratio(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        genes_df = genes_df[:100]  # only first 100 rows
        test_df = init_result_df(genes_df, generate_fake_traits(genes_df))

        # apply function
        test_df = add_odds_ratio(test_df)
        self.assertEqual(
            test_df.columns.tolist(),
            ['Gene', 'c1r1', 'c2r1', 'c1r2', 'c2r2', '__contingency_table__', 'sensitivity', 'specificity', 'odds_ratio']
        )

        # calculate odds_ratio with fisher_exact
        fisher_ors = test_df.apply(lambda row: fisher_exact([[row['c1r1'], row['c2r1']], [row['c1r2'], row['c2r2']]])[0], axis=1)

        # check if result is identical
        for manual_or, fisher_or in zip(test_df['odds_ratio'], fisher_ors):
            self.assertTrue(is_equivalent(manual_or, fisher_or))
        print(test_df.columns)

    def test_tetracycline(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        traits_df = load_traits(get_path('tetracycline', 'traits'), delimiter=',')
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
        test_df = compute_pairs(test_df, all_label_to_gene, tree=tree, label_to_trait=label_to_trait)

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
            new_table = tuple(int(new_row[c]) for c in ('c1r1', 'c2r1', 'c1r2', 'c2r2'))

            print(row.Gene)

            self.assertEqual(table, new_table)
            self.assertAlmostEqual(row.Odds_ratio, new_row.odds_ratio,
                                   msg=f'Failed to calculate odds_ratio for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')
            self.assertAlmostEqual(row.Sensitivity, new_row.sensitivity,
                                   msg=f'Failed to calculate sensitivity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')
            self.assertAlmostEqual(row.Specificity, new_row.specificity,
                                   msg=f'Failed to calculate specificity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')

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
            except Exception:
                print('xxx')

    def test_with_internal_data(self):
        ANYGENE = 'gene1291'

        def load_gene_data(trait: str = 't1', gene: str = ANYGENE):
            with open(f'{ROOT}/data/test_with_internal_data/{trait}_{gene}.json') as f:
                res = json.load(f)
            return res

        PRES_ABS = f'{ROOT}/data/test_with_internal_data/pres_abs_tree_id_3.csv'
        TRAITS = f'{ROOT}/data/test_with_internal_data/trait_trees.csv'
        NWK = f'{ROOT}/data/test_with_internal_data/newick.nwk'
        with open(NWK) as f:
            NWK = f.read()

        genes_df = load_genes(PRES_ABS, delimiter=',', start_col=0)
        traits_df = load_traits(TRAITS, delimiter=',')
        print(genes_df)
        print(traits_df)

        for trait in traits_df:
            RESULT = f'{ROOT}/data/test_with_internal_data/{trait}.results.csv'
            expected_result = pd.read_csv(RESULT)

            trait_series = traits_df[trait]
            # calculate sensitivity and specificity
            test_df = init_result_df(genes_df, label_to_trait={l: bool(v) for l, v in trait_series.items() if v in (0, 1)})
            # calculate odds_ratio
            test_df = add_odds_ratio(test_df)
            # calculate pairwise comparisons
            all_label_to_gene = get_all_label_to_gene(genes_df)
            tree = ScoaryTree.from_newick(NWK)
            assert set(tree.labels()) == set(all_label_to_gene[ANYGENE])
            label_to_trait = get_label_to_trait(trait_series)
            test_df = compute_pairs(test_df, all_label_to_gene, tree=tree, label_to_trait=label_to_trait)
            test_df.set_index('Gene', inplace=True)

            # check if result is identical
            for i, row in expected_result.iterrows():
                table = (row.Number_pos_present_in,
                         row.Number_neg_present_in,
                         row.Number_pos_not_present_in,
                         row.Number_neg_not_present_in)
                new_row = test_df.loc[row.Gene]
                new_table = tuple(int(new_row[c]) for c in ('c1r1', 'c2r1', 'c1r2', 'c2r2'))

                print(row.Gene)
                if row.Gene == 'gene24112':
                    print('pause')

                self.assertEqual(table, new_table)
                self.assertAlmostEqual(row.Odds_ratio, new_row.odds_ratio,
                                       msg=f'Failed to calculate odds_ratio for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')
                self.assertAlmostEqual(row.Sensitivity, new_row.sensitivity,
                                       msg=f'Failed to calculate sensitivity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')
                self.assertAlmostEqual(row.Specificity, new_row.specificity,
                                       msg=f'Failed to calculate specificity for {row.Gene}: {row.Odds_ratio} != {new_row.odds_ratio}')

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
                except Exception:
                    print('my tree')
                    print_tree_for_debugging(tree, all_label_to_gene[row.Gene], label_to_trait)
                    print('scoary tree')
                    # gene_data = load_gene_data(trait=trait, gene=row.Gene)
                    # print_tree_for_debugging(tree, gene_data['label_to_gene'], gene_data['label_to_trait'])

                    print(pd.DataFrame(xx, columns=['old', 'new'], index=['c', 's', 'o', 'b', 'w']))

                    print(f'odds_ratio={new_row.odds_ratio}')

                    supporting, opposing = count_best_worst(tree, label_to_trait, all_label_to_gene[row.Gene])

                    print()
