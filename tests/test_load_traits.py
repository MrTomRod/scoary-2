from init_tests import *
from scoary.load_traits import load_numeric, load_binary, apply_kmeans, apply_gm, binarize, load_traits

traits_bin = get_path('tetracycline', 'traits')
traits_num = get_path('tetracycline', 'traits-numeric')


class Test(TestCase):
    def test_load_binary(self):
        binary_df = load_binary(traits=traits_bin, delimiter=',')
        print(binary_df)

    def test_load_numeric(self):
        numeric_df = load_numeric(traits=traits_num, delimiter=',')
        print(numeric_df)

    def test_binarize_kmeans(self):
        numeric_df = load_numeric(traits=traits_num, delimiter=',')
        for alternative in ['skip', 'kmeans']:
            for cutoff in [.5, .7, .9]:
                for covar_type in ['full', 'tied', 'diag', 'spherical']:
                    binary_df = binarize(
                        numeric_df, method='kmeans', random_state=42, n_cpus=1,
                        cutoff=cutoff, covariance_type=covar_type,
                        alternative=alternative, outdir=None
                    )

    def test_binarize_gaussian_nonconverging(self):
        # Tetracycline trait cannot be binarized with cutoff=0.999
        numeric_df = load_numeric(traits=traits_num, delimiter=',')
        for method, n_expected_columns in [('gaussian', 1), ('kmeans', 2)]:
            binary_df = binarize(
                numeric_df, method=method, random_state=42, n_cpus=1,
                cutoff=0.9998, covariance_type='full',
                alternative='skip', outdir=None
            )
            self.assertEqual(n_expected_columns, len(binary_df.columns),
                             f'{method=}; {n_expected_columns=}; {binary_df.columns=}')

    def test_illegal(self):
        with self.assertRaises(AssertionError):
            numeric_df, traits_df = load_traits(traits_num, trait_data_type=f'gaussian:0.4999', random_state=42)
        with self.assertRaises(AssertionError):
            numeric_df, traits_df = load_traits(traits_num, trait_data_type=f'gaussian:1', random_state=42)
        with self.assertRaises(AssertionError):
            # fails because no traits can be binarized. Certainty is never high enough.
            numeric_df, traits_df = load_traits(traits_num, trait_data_type=f'gaussian:.999999999999', random_state=42)

    def test_multiprocessing(self):
        for n_cpus in [1, 5]:
            numeric_df, traits_df = load_traits(
                get_path('new_ds', 'traits-lc'),
                trait_data_type='gaussian:skip:\t',
                ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture',
                n_cpus=n_cpus, limit_traits=(0, 10),
                outdir=f'{ROOT}/TEST_OUTPUT'
            )
