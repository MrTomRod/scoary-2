from os.path import dirname
import json
from unittest import TestCase
from scipy.spatial import distance
from scipy.stats import fisher_exact

# set up logging
import logging

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

from scoary.scoary import *

DATA = {
    'generated': {
        'genes': 'Gene_presence_absence.csv',
        'traits': 'Trait.csv'
    },
    'tetracycline': {
        'genes': 'Gene_presence_absence.csv',
        'traits': 'Tetracycline_resistance.csv',
        'restrict_to': 'Restrict_to.csv',
        'tree': 'ExampleTree.nwk',
        'treelist': 'expected_result.json',
        'tdm': 'tetracycline_TDM.csv'
    }
}
VCF_data = {
    'vcf': {
        'vcf': 'Example.vcf',
        'traits': 'ExampleVCFTrait.csv'
    }
}

ROOT = dirname(dirname(__file__))


def get_path(ds: str, key: str, data=DATA):
    return f'{ROOT}/data/{ds}/{data[ds][key]}'


def is_equivalent_tree(a, b) -> bool:
    if type(a) is str or type(b) is str:
        return a == b
    else:
        return (is_equivalent_tree(a[0], b[0]) and is_equivalent_tree(a[1], b[1])) or (
                is_equivalent_tree(a[0], b[1]) and is_equivalent_tree(a[1], b[0]))


def is_equivalent(a, b):
    if np.isinf(a) and np.isinf(b):
        return True
    if np.isnan(a) and np.isnan(b):
        return True
    return np.isclose(a, b)


class Test(TestCase):
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

    def test_create_tree(self):
        """
        Check if old scoary generates the equivalent tree based on genes presence/absence
        """
        with open(get_path('tetracycline', 'treelist')) as f:
            expected_result = json.load(f)['as_list']

        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        res = tree_from_presence_absence(genes_df)

        self.assertTrue(is_equivalent_tree(res, expected_result))

    def test_load_tree(self):
        with open(get_path('tetracycline', 'treelist')) as f:
            newick = json.load(f)['as_newick']

        tree_1 = tree_from_file(newick=newick)
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        tree_2 = tree_from_presence_absence(genes_df)

        self.assertTrue(is_equivalent_tree(tree_1, tree_2))

    def test_create_test_df(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        test_df = create_test_df(genes_df, trait_pos=set(genes_df.columns[:10]), trait_neg=set(genes_df.columns[90:]))
        self.assertEqual(test_df.columns.tolist(), ['c1r1', 'c2r1', 'c1r2', 'c2r2'])

    def test_contingency_test(self):
        for CachedContingency in CachedFisher, CachedBoschloo:
            cc = CachedContingency()
            genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
            test_df = create_test_df(genes_df, trait_pos=set(genes_df.columns[:11]), trait_neg=set(genes_df.columns[89:]))
            test_df_2 = perform_contingency_test(test_df=test_df, cc=cc)
            self.assertEqual(test_df_2.columns.tolist(), ['c1r1', 'c2r1', 'c1r2', 'c2r2', f'pval_{cc.function_name}'])
            print(f"Done: {cc} -> minpval={test_df_2[f'pval_{cc.function_name}'].min()}")

    def test_odds_ratio(self):
        genes_df = load_genes(get_path('tetracycline', 'genes'), delimiter=',', start_col=13)
        genes_df = genes_df[:100]  # only first 100 rows
        test_df = create_test_df(genes_df, trait_pos=set(genes_df.columns[:11]), trait_neg=set(genes_df.columns[89:]))

        # apply function
        test_df = add_odds_ratio(test_df)
        self.assertEqual(test_df.columns.tolist(), ['c1r1', 'c2r1', 'c1r2', 'c2r2', 'odds_ratio'])

        # calculate odds_ratio with fisher_exact
        fisher_ors = test_df.apply(lambda row: fisher_exact([[row['c1r1'], row['c2r1']], [row['c1r2'], row['c2r2']]])[0], axis=1)

        # check if result is identical
        for manual_or, fisher_or in zip(test_df['odds_ratio'], fisher_ors):
            self.assertTrue(is_equivalent(manual_or, fisher_or))
        print(test_df.columns)
