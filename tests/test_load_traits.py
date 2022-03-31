from init_tests import *
from scoary.load_traits import binarize

traits_file = get_path('tetracycline', 'traits-numeric')


class Test(TestCase):
    def test_binarize_kmeans(self):
        res = binarize(traits_file, trait_data_type='kmeans')

    def test_binarize_gaussian(self):
        for alternative in ['skip', 'kmeans']:
            for cutoff in [.5, .7, .9]:
                for covar_type in ['full', 'tied', 'diag', 'spherical']:
                    numeric_df, traits_df = binarize(
                        traits_file, trait_data_type=f'gaussian:{alternative}:{cutoff}:{covar_type}', random_state=42
                    )

    def test_binarize_gaussian_nonconverging(self):
        # Tetracycline trait cannot be binarized with cutoff=0.999
        numeric_df, traits_df = binarize(traits_file, trait_data_type=f'gaussian:skip:0.999', random_state=42)
        print(traits_df)
        self.assertEqual(len(traits_df.columns), 1)

        # use kmeans if gaussian doesn't work
        numeric_df, traits_df = binarize(traits_file, trait_data_type=f'gaussian:kmeans:0.999', random_state=42)
        self.assertEqual(len(traits_df.columns), 2)

    def test_illegal(self):
        with self.assertRaises(AssertionError):
            numeric_df, traits_df = binarize(traits_file, trait_data_type=f'gaussian:0.4999', random_state=42)
        with self.assertRaises(AssertionError):
            numeric_df, traits_df = binarize(traits_file, trait_data_type=f'gaussian:1', random_state=42)
        with self.assertRaises(AssertionError):
            # fails because no traits can be binarized. Certainty is never high enough.
            numeric_df, traits_df = binarize(traits_file, trait_data_type=f'gaussian:.999999999999', random_state=42)

    def test_new_ds(self):
        numeric_df, traits_df = binarize(
            get_path('new_ds', 'traits-lc'),
            trait_data_type=f'gaussian:skip:\t',
            ignore='Starter-only-5A,FAMIX,Starter-only-10,Starter-only-7,mixture'
        )
        print(numeric_df)
        print(traits_df)
