import shutil

from init_tests import *

from scoary.scoary import *


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
            n_cpus=1,
            outdir=self.tempdir
        )

    def test_scoary_multi_threaded(self):
        scoary(
            genes=get_path('tetracycline', 'genes'),
            traits=get_path('tetracycline', 'traits'),
            n_permut=1000,
            n_cpus=4,
            outdir=self.tempdir
        )

    def test_scoary_long(self):
        scoary(
            genes=get_path('new_ds', 'genes-hog'),
            gene_data_type='gene-list:\t',
            traits=get_path('new_ds', 'traits-lc-binary'),
            trait_data_type='binary:\t',
            n_permut=200,
            # ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
            restrict_to='FAM1414-i1-1.1,FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            random_state=42,
            n_cpus=7,
            outdir=self.tempdir,
            limit_traits=(0, 20),
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
            restrict_to='FAM1414-i1-1.1,FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,'
                        'FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,'
                        'FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,'
                        'FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,'
                        'FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,'
                        'FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,'
                        'FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,'
                        'FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,'
                        'FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            random_state=42,
            n_cpus=7,
            outdir=self.tempdir,
            limit_traits=(0, 20),
            pairwise=False
        )

    def test_scoary_real(self):
        scoary(
            # multiple_testing_fisher='bonferroni:0.01',
            genes=get_path('full_ds', 'genes'),
            gene_info=get_path('full_ds', 'gene-info'),
            gene_data_type='gene-list:\t',
            traits=get_path('full_ds', 'traits'),
            trait_data_type=f'gaussian:skip:\t:tied',  # {'tied', 'full', 'diag', 'spherical'}
            trait_info=get_path('full_ds', 'trait-info'),
            isolate_info=get_path('full_ds', 'isolate-info'),
            n_permut=300,
            random_state=42,
            n_cpus=8,
            restrict_to='FAM14177-p1-1.1,FAM14184-i1-1.1,FAM14193-i1-1.1,FAM14197-i1-1.1,FAM14217-p1-1.1,FAM14221-p1-1.1,FAM14222-p1-1.1,FAM1414-i1-1.1,FAM15061-i1-1.1,FAM15078-i1-1.1,FAM15113-i1-1.1,FAM15170-i1-1.1,FAM15190-i1-1.1,FAM15192-i1-1.1,FAM15300-i1-1.1,FAM15333-i1-1.1,FAM15346-i1-1.1,FAM15347-i1-1.1,FAM15381-i1-1.1,FAM15407-i1-1.1,FAM19015-i1-1.1,FAM19016-i1-1.1,FAM19020-i1-1.1,FAM19022-i1-1.1,FAM19023-i1-1.1,FAM19024-p1-1.1,FAM19025-p1-1.1,FAM19030-i2-1.1,FAM19031-i2-1.1,FAM19034-i1-1.1,FAM22019-i1-1.1,FAM22020-i1-1.1,FAM22021-p1-1.1,FAM23848-i1-1.1,FAM23852-i1-1.1,FAM23853-i1-1.1,FAM23855-i1-1.1,FAM23864-i1-1.1,FAM23867-i1-1.1,FAM23868-i1-1.1,FAM23869-i1-1.1,FAM23870-i1-1.1,FAM23877-p1-1.1,FAM24252-i1-1.1',
            limit_traits=(2330, 2340),
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
        )

    def test_same_hemming_result(self):
        """
        Check if old scoary generates the same data (hamming similarity matrix)
        """
        _, genes_df = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count', ignore=roary_ignore)
        tdm_new = pd.DataFrame(distance.squareform(distance.pdist(genes_df.T, 'hamming')))
        tdm_old = np.flip(pd.read_csv(get_path('tetracycline', 'tdm'), index_col=0).values)  # has to be flipped
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
