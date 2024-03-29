import os
from subprocess import call
from init_tests import *
from scoary.final_overview import create_final_overview
from scoary.load_traits import load_binary
from scoary.utils import pd, AnalyzeTraitNamespace

REPLACE_COPIES_WITH_SYMLINKS = True


def replace_copies_with_symlinks():
    def repl(fn, relpath='../..', subdir='app'):
        src = f'{relpath}/scoary/templates/{fn}'
        target = f'../TEST_OUTPUT/{subdir}/{fn}'
        if os.path.isfile(target):
            os.remove(target)
        os.symlink(src=src, dst=target)

    for file in ['trait.html', 'overview.html']:
        repl(file, relpath='..', subdir='')
    for file in ['config.json', 'favicon.svg', 'overview.css', 'overview.js', 'trait.css', 'trait.js']:
        repl(file)


class Test(TestCase):
    def setUp(self) -> None:
        self.temp_dir = get_tempdir_path()
        self.fake_ns = AnalyzeTraitNamespace()
        self.fake_ns.outdir = self.temp_dir
        self.fake_ns.trait_info_df = None

        os.makedirs(self.temp_dir, exist_ok=True)
        call(f'rm -rf {self.temp_dir}/*', shell=True)
        for dir_ in ['app', 'traits', 'logs']:
            os.makedirs(f'{self.temp_dir}/{dir_}', exist_ok=True)

    def tearDown(self) -> None:
        if REPLACE_COPIES_WITH_SYMLINKS:
            replace_copies_with_symlinks()
        print(f'Open file://{self.temp_dir} to see the result!')
        print(f'To clean up, run "rm -r {self.temp_dir}"')

    def test_simple(self):
        summary_df = pd.DataFrame(**{'index': ['Compound_242', 'Compound_267', 'Compound_286'],
                                     'columns': ['best_fisher_p', 'best_fisher_q', 'best_empirical_p', 'best_fq*ep'],
                                     'data': [[0.574065934065931, 0.438405797101457, 0.03596403596403, 1.576684e-02],
                                              [0.432940190858691, 0.266793137470672, 0.13386613386613, 3.571457e-02],
                                              [0.194418465932588, 7.98120572982e-08, 0.02097902097902, 1.674379e-09]]})
        traits_df = pd.DataFrame(**{
            'index': ['FAM10789-i1-1.1', 'FAM1079-i1-1.1', 'FAM10792-i1-1.1', 'FAM11142-i1-1.1', 'FAM11194-i1-1.1',
                      'FAM11199-i1-1.1', 'FAM11206-i1-1.1', 'FAM1233-i1-1.1', 'FAM1301-i1-1.1', 'FAM13493-i1-1.1'],
            'columns': ['Compound_242', 'Compound_267', 'Compound_286'],
            'data': [[pd.NA, True, True], [pd.NA, pd.NA, True], [pd.NA, pd.NA, True], [pd.NA, False, False],
                     [pd.NA, True, True], [pd.NA, pd.NA, True], [True, pd.NA, True], [pd.NA, False, True],
                     [False, pd.NA, True], [pd.NA, pd.NA, True]]},
                                 dtype='boolean')
        self.fake_ns.traits_df = traits_df  # load_binary('../data/new_ds/LC-binary.tsv', '\t')
        create_final_overview(summary_df=summary_df, ns=self.fake_ns)

    def test_larger(self):
        self.fake_ns.traits_df = load_binary('../data/new_ds/LC-binary.tsv', '\t')
        summary_df = pd.DataFrame(index=self.fake_ns.traits_df.columns)
        for col in ['best_fisher_p', 'best_fisher_q', 'best_empirical_p']:
            summary_df[col] = np.random.rand(1, len(self.fake_ns.traits_df.columns))[0]
        create_final_overview(summary_df=summary_df, ns=self.fake_ns)

    def test_largest(self):
        # This function was used to determine the desired recursion limit in plot_dendrogram
        n_traits, n_isolates = 100, 44  # 10000, 44
        self.fake_ns.traits_df = pd.DataFrame(
            np.random.rand(n_isolates, n_traits) > 0.5,
            index=[f'I{i}' for i in range(n_isolates)],
            columns=[f'T{i}' for i in range(n_traits)],
        )
        summary_df = pd.DataFrame(index=self.fake_ns.traits_df.columns)
        for col in ['best_fisher_p', 'best_fisher_q', 'best_empirical_p']:
            summary_df[col] = np.random.rand(1, len(self.fake_ns.traits_df.columns))[0]

        create_final_overview(summary_df=summary_df, ns=self.fake_ns)

    def test_real(self):
        # Find the best way of plotting dendrogram
        # summary_df = pd.read_csv(f'../data/summary__.tsv', sep='\t', index_col=0)
        summary_df = pd.read_csv(f'../TMP/TEST_OUTPUT_real_restricted/summary.tsv', sep='\t', index_col=0)
        self.fake_ns.traits_df = pd.read_csv('../data/traits_df.tsv', sep='\t', index_col=0, dtype='str') == 'True'
        create_final_overview(
            summary_df=summary_df,
            traits_df=self.fake_ns.traits_df,
            outdir=self.fake_ns.outdir
        )


    def test_understand_jaccard(self):
        from scipy.spatial.distance import cdist

        a = np.array([[0, 0, 0, 0, 0, 0, 1, 0, -1, 0],
                      [1, 0, 0, -1, 1, 0, 0, -1, 0, 0],
                      [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, ],
                      [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, ],
                      [1, 1, 1, -1, 1, 1, 1, 1, 1, 1]], dtype=int)
        print('input:\n', a)

        def conf(arr):
            aa, nn = arr.shape
            res = np.zeros(shape=(aa, aa))

            for a in range(aa):
                for b in range(aa):
                    if a == b:
                        break

                    n = nn - sum(np.logical_and(arr[a] == 0, arr[b] == 0))

                    # x = np.abs(arr[a] - arr[b]) / 2
                    # y = np.abs(arr[a] - (0 - arr[b])) / 2
                    # r = min(x.sum(), y.sum()) / n
                    # print(a, b, x, y, n, r)

                    x = arr[a] != arr[b]
                    y = arr[a] != (0 - arr[b])
                    r = min(x.sum(), y.sum()) / n
                    print(a, b, x, y, n, r)

                    res[a, b] = r
                    res[b, a] = r
            print('res')
            print(res)

        conf(a)

        def x(a):
            a = np.nan_to_num(a, nan=0.5)
            b = 0 - a

            d1 = cdist(a, a, metric='jaccard')
            d2 = cdist(a, b, metric='jaccard')
            d = np.minimum(d1, d2)

            # print('a')
            # print(a)
            # print('b')
            # print(b)
            # print('d1')
            # print(d1)
            # print('d2')
            # print(d2)
            print('d')
            print(d)

        x(a)


class ReplaceCopiesWithSymlinks(TestCase):
    def test_replace_copies_with_symlinks(self):
        replace_copies_with_symlinks()
