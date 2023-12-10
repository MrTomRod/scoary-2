import shutil

from init_tests import *

from scoary.scoary import *

os.environ['MGWAS_LINK_ONLY'] = 'true'

RESTRICT_TO = 'FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,' \
              'FAM14222-p1-1.1,FAM1414-i1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,' \
              'FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,' \
              'FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,' \
              'FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,' \
              'FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,' \
              'FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,' \
              'FAM23877-p1-1.1,FAM24252-i1-1.1'


class TestScoary(TestCase):
    def setUp(self) -> None:
        self.tempdir = get_tempdir_path()
        if os.path.isdir(self.tempdir):
            shutil.rmtree(self.tempdir)

    def test_scoary_single_threaded(self):
        scoary(
            trait_wise_correction=False,
            genes='../data/tetracycline/Gene_presence_absence.csv',
            traits='../data/tetracycline/Tetracycline_resistance.csv',
            n_permut=1000,
            multiple_testing='fdr_bh:0.5',
            n_cpus=1,
            outdir=self.tempdir
        )

    def test_scoary_multi_threaded(self):
        scoary(
            trait_wise_correction=True,
            genes='../data/tetracycline/Gene_presence_absence.csv',
            traits='../data/tetracycline/Tetracycline_resistance.csv',
            n_permut=200,
            n_cpus=4,
            outdir=self.tempdir,
            multiple_testing='native:0.05'
        )

    def test_scoary_gene_info(self):
        scoary(
            genes='../data/tetracycline/Gene_presence_absence.csv',
            gene_info='../data/tetracycline/gene-info.tsv',
            traits='../data/tetracycline/Tetracycline_resistance.csv',
            n_permut=10000,
            n_cpus=1,
            outdir=self.tempdir
        )

    def test_scoary_long_binary(self):
        scoary(
            trait_wise_correction=True,
            multiple_testing='fdr_bh:0.6',
            linkage_method='average',
            genes='../data/new_ds/N0.tsv',
            gene_data_type='gene-list:\t',
            traits='../data/new_ds/LC-binary.tsv',
            trait_data_type='binary:\t',
            n_permut=200,
            # ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
            restrict_to=RESTRICT_TO,
            random_state=42,
            n_cpus=7,
            outdir=self.tempdir,
            # limit_traits=(0, 20),
            limit_traits=(320, 340),
            max_genes=100
        )

    def test_scoary_long_numeric(self):
        scoary(
            multiple_testing='fdr_bh:0.3',
            genes='../data/new_ds/N0.tsv',
            gene_info='../data/new_ds/N0_best_names.tsv',
            gene_data_type='gene-list:\t',
            traits='../data/new_ds/LC.tsv',
            trait_data_type='gaussian:skip:\t:tied',
            trait_info='../data/new_ds/LC-meta.tsv',
            isolate_info='../data/new_ds/isolate-meta.tsv',
            n_permut=200,
            # ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
            restrict_to='FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM1414-i1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            random_state=42,
            n_cpus=7,
            outdir=self.tempdir,
            limit_traits=(0, 200),
            pairwise=True
        )

    def test_scoary_gauss_kmeans(self):
        scoary(
            genes='../data/new_ds/N0.tsv',
            gene_info='../data/new_ds/N0_best_names.tsv',
            gene_data_type='gene-list:\t',
            traits='../data/new_ds/LC.tsv',
            trait_data_type=f'gaussian:kmeans:\t',
            trait_info='../data/new_ds/LC-meta.tsv',
            isolate_info='../data/new_ds/isolate-meta.tsv',
            n_permut=200,
            restrict_to=RESTRICT_TO,
            random_state=42,
            n_cpus=7,
            outdir=self.tempdir,
            # limit_traits=(0, 100),
            # pairwise=False
        )

    def test_scoary_full(self):
        scoary(
            multiple_testing='bonferroni:0.1',
            genes='../data/full_ds/N0.tsv',
            gene_info='../data/full_ds/N0_best_names.tsv',
            gene_data_type='gene-list:\t',
            traits='../data/full_ds/traits.tsv',
            trait_data_type=f'gaussian:skip:\t:tied',  # {'tied', 'full', 'diag', 'spherical'}
            trait_info='../data/full_ds/trait_info.tsv',
            isolate_info='../data/full_ds/isolate_info.tsv',
            n_permut=600,
            random_state=42,
            n_cpus=8,
            n_cpus_binarization=1,
            restrict_to=RESTRICT_TO,
            max_genes=50,
            # limit_traits=(12377, 12378),
            limit_traits=(2000, 2400),
            # limit_traits=(2330, 2340),
            worst_cutoff=0.1,
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
            n_cpus=1,
            outdir=self.tempdir,
            multiple_testing='native:0.05',
        )

    def test_scoary_jacordova(self):
        scoary(
            genes='../data/jacordova/GeneCount_Scoary_Ecoli.txt',
            gene_data_type='gene-count:\t',
            traits='../data/jacordova/Ecoli_traits.txt',
            trait_data_type='gaussian:kmeans:\t',  # {'tied', 'full', 'diag', 'spherical'}
            n_permut=1000,
            random_state=42,
            n_cpus=1,
            outdir=self.tempdir,
            multiple_testing='native:0.05',
        )

    def test_same_hemming_result(self):
        """
        Check if old scoary generates the same data (hamming similarity matrix)
        """
        _, genes_df = load_genes('../data/tetracycline/Gene_presence_absence.csv', gene_data_type='gene-count', ignore=roary_ignore)
        tdm_new = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')))
        tdm_old = np.flip(pd.read_csv('../data/tetracycline/tetracycline_TDM.csv', index_col=0).values)  # has to be flipped
        np.fill_diagonal(tdm_old, 0)  # diagonal should be 0, not 1
        tdm_old = pd.DataFrame(tdm_old)
        self.assertTrue(np.isclose(tdm_old, tdm_new).all())

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
            n_cpus=4,
            outdir=self.tempdir
        )

    def test_scoary_roary_gene_list(self):
        # GitHub issue #5
        # scoary(
        #     genes=get_path('roary-list', 'genes'),
        #     traits=get_path('roary-list', 'traits'),
        #     gene_data_type='gene-list:,',
        #     n_permut=1000,
        #     multiple_testing='native:0.05',
        #     n_cpus=1,
        #     outdir=self.tempdir
        # )
        scoary(
            genes='../data/roary-list/gene_presence_absence-b.csv',
            traits='../data/roary-list/traits-b.csv',
            gene_data_type='gene-list:,',
            n_permut=1000,
            multiple_testing='native:0.05',
            n_cpus=1,
            outdir=self.tempdir
        )

    def test_scoary_pyseer(self):
        scoary(
            genes='../data/pyseer/gene_presence_absence.Rtab',
            traits='../data/pyseer/resistances.pheno',
            gene_data_type='gene-count:\t',
            trait_data_type='binary:\t',
            multiple_testing='bonferroni:0.05',
            n_cpus=1,
            outdir=self.tempdir,
            pairwise=False
        )
