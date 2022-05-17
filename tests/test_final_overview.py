import os
from subprocess import call
from init_tests import *
from scoary.final_overview import create_final_overview
from scoary.load_traits import load_binary
from scoary.utils import pd, AnalyzeTraitNamespace


class Test(TestCase):
    def setUp(self) -> None:
        self.temp_dir = get_tempdir_path()
        self.fake_ns = AnalyzeTraitNamespace()
        self.fake_ns.outdir = self.temp_dir
        self.fake_ns.trait_info_df = None

        os.makedirs(self.temp_dir, exist_ok=True)
        call(f'rm -rf {self.temp_dir}/*', shell=True)

    def tearDown(self) -> None:
        print(f'Open file://{self.temp_dir} to see the result!')
        print(f'To clean up, run "rm -r {self.temp_dir}"')

    def test_simple(self):
        summary_df = pd.DataFrame(**{'index': ['Compound_242', 'Compound_267', 'Compound_286'],
                                      'columns': ['min_pval', 'min_qval', 'min_pval_empirical', 'min_qval_empirical'],
                                      'data': [[0.574065934065931, 0.43840579710145217, 0.03596403596403597,
                                                0.03596403596403597],
                                               [0.4329401908586914, 0.2667931374706272, 0.13386613386613386,
                                                0.13386613386613386],
                                               [0.19441846593258844, 7.981205729820185e-08, 0.02097902097902098,
                                                0.6913086913086913]]})
        traits_df = pd.DataFrame(**{
            'index': ['FAM10789-i1-1.1', 'FAM1079-i1-1.1', 'FAM10792-i1-1.1', 'FAM11142-i1-1.1', 'FAM11194-i1-1.1',
                      'FAM11199-i1-1.1', 'FAM11206-i1-1.1', 'FAM1233-i1-1.1', 'FAM1301-i1-1.1', 'FAM13493-i1-1.1'],
            'columns': ['Compound_242', 'Compound_267', 'Compound_286'],
            'data': [[pd.NA, True, True], [pd.NA, pd.NA, True], [pd.NA, pd.NA, True], [pd.NA, False, False],
                     [pd.NA, True, True], [pd.NA, pd.NA, True], [True, pd.NA, True], [pd.NA, False, True],
                     [False, pd.NA, True], [pd.NA, pd.NA, True]]},
                                 dtype='boolean')
        self.fake_ns.traits_df = traits_df  # load_binary(get_path('new_ds', 'traits-lc-binary'), '\t')
        create_final_overview(summary_df=summary_df, ns=self.fake_ns)

    def test_larger(self):
        self.fake_ns.traits_df = load_binary(get_path('new_ds', 'traits-lc-binary'), '\t')
        summary_df = pd.DataFrame(index=self.fake_ns.traits_df.columns)
        for col in ['min_pval', 'min_qval', 'min_pval_empirical', 'min_qval_empirical']:
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
        for col in ['min_pval', 'min_qval', 'min_pval_empirical', 'min_qval_empirical']:
            summary_df[col] = np.random.rand(1, len(self.fake_ns.traits_df.columns))[0]

        create_final_overview(summary_df=summary_df, ns=self.fake_ns)

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
