import shutil

from init_tests import *

from scoary.scoary import *


def generate_fake_traits(genes_df: pd.DataFrame) -> {str: bool}:
    label_to_trait = {}
    label_to_trait.update({l: True for l in genes_df.columns[:11]})
    label_to_trait.update({l: False for l in genes_df.columns[89:]})
    return label_to_trait


class TestScoary(TestCase):
    def setUp(self) -> None:
        self.tempdir = get_tempdir_path()
        if os.path.isdir(self.tempdir):
            shutil.rmtree(self.tempdir)

    def test_scoary_single_threaded(self):
        scoary(
            genes=get_path('tetracycline', 'genes'),
            traits=get_path('tetracycline', 'traits'),
            n_permut=1000,
            threads=1,
            outdir=self.tempdir
        )

    def test_scoary_multi_threaded(self):
        scoary(
            genes=get_path('tetracycline', 'genes'),
            traits=get_path('tetracycline', 'traits'),
            n_permut=1000,
            threads=4,
            outdir=self.tempdir
        )

    def test_scoary_long(self):
        scoary(
            genes=get_path('new_ds', 'genes-hog'),
            gene_data_type='gene-list:2',
            traits=get_path('new_ds', 'traits-lc-binary'),
            trait_data_type='binary:\t',
            n_permut=200,
            # ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
            restrict_to='FAM1414-i1-1.1,FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            random_state=42,
            threads=7,
            outdir=self.tempdir,
            limit_to_n_traits=10,
        )

    def test_scoary_long_numeric(self):
        scoary(
            genes=get_path('new_ds', 'genes-hog'),
            gene_info=get_path('new_ds', 'genes-hog-info'),
            gene_data_type='gene-list:2',
            traits=get_path('new_ds', 'traits-lc'),
            trait_data_type=f'gaussian:skip:\t:tied',  # {'tied', 'full', 'diag', 'spherical'}
            trait_info=get_path('new_ds', 'traits-lc-meta'),
            isolate_info=get_path('new_ds', 'isolate-meta'),
            n_permut=200,
            # ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
            restrict_to='FAM1414-i1-1.1,FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            random_state=42,
            threads=7,
            outdir=self.tempdir,
            # no_pairwise=True
        )

    def test_scoary_real(self):
        scoary(
            genes=get_path('full_ds', 'genes'),
            gene_info=get_path('full_ds', 'gene-info'),
            gene_data_type='gene-list:\t',
            traits=get_path('full_ds', 'traits'),
            trait_data_type=f'gaussian:skip:\t:tied',  # {'tied', 'full', 'diag', 'spherical'}
            trait_info=get_path('full_ds', 'trait-info'),
            isolate_info=get_path('full_ds', 'isolate-info'),
            n_permut=300,
            random_state=42,
            threads=7,
            restrict_to='FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM1414-i1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            limit_traits=(0, 10),
            outdir=self.tempdir,
        )

    def test_scoary_marco(self):
        scoary(
            genes='../data/marco/Orthogroups.tsv',
            gene_data_type='gene-list:\t',
            traits='../data/marco/traits.tsv',
            trait_data_type='binary: ',  # {'tied', 'full', 'diag', 'spherical'}
            n_permut=1000,
            random_state=42,
            threads=1,
            outdir=self.tempdir,
        )

    def test_same_hemming_result(self):
        """
        Check if old scoary generates the same data (hamming similarity matrix)
        """
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
        tdm_new = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')))
        tdm_old = np.flip(pd.read_csv(get_path('tetracycline', 'tdm'), index_col=0).values)  # has to be flipped
        np.fill_diagonal(tdm_old, 0)  # diagonal should be 0, not 1
        tdm_old = pd.DataFrame(tdm_old)
        self.assertTrue(np.isclose(tdm_old, tdm_new).all())

    def test_create_result_df(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        self.assertEqual(
            result_df.columns.tolist(),
            ['Gene', 'g+t+', 'g+t-', 'g-t+', 'g-t-', '__contingency_table__', 'sensitivity', 'specificity']
        )

    def test_contingency_test(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
        result_df = init_result_df(genes_df, label_to_trait=generate_fake_traits(genes_df))
        test_df = create_test_df(result_df=result_df)
        self.assertEqual(['__contingency_table__', 'pval'], test_df.columns.tolist())
        print(f"Done: minpval={test_df.pval.min()}")

    def test_odds_ratio(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
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
        print(test_df.columns)

    def test_tetracycline(self):
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=tetr_ignore)
        _, traits_df = load_traits(get_path('tetracycline', 'traits'), delimiter=',')
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
            except Exception as e:
                print(i, row.Gene, xx)
                self.fail(msg=str(e))

    def test_recursion_depth(self):
        strains = [f'strain_{i}' for i in range(13000)]
        genes = [f'gene_{i}' for i in range(100)]
        traits = [f'trait_{i}' for i in range(4)]
        genes_df = pd.DataFrame(
            np.random.randint(
                low=0, high=2, size=(len(genes), len(strains))
            ), index=genes, columns=strains
        )
        traits_df = pd.DataFrame(
            np.random.randint(
                low=0, high=2, size=(len(strains), len(traits))
            ), index=strains, columns=traits
        )
        genes_df.to_csv('../data/huge_ds/genes.tsv', sep='\t')
        traits_df.to_csv('../data/huge_ds/traits.tsv', sep='\t')
        # Calculating tree is very slow, but it works.
        with open('../data/huge_ds/tree.nwk', 'w') as f:
            f.write('(' * (len(strains) - 1))
            f.write(strains[0])
            f.write(',')
            f.write('),'.join(strains[1:]))
            f.write(');')

        scoary(
            genes='../data/huge_ds/genes.tsv',
            traits='../data/huge_ds/traits.tsv',
            trait_data_type='binary:\t',
            gene_data_type='gene-count:\t',
            newicktree='../data/huge_ds/tree.nwk',
            n_permut=1000,
            threads=4,
            outdir=self.tempdir
        )
