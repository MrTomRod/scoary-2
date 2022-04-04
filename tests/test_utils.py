from init_tests import *
from scoary.utils import *

logger = logging.getLogger('TEST_LOGGER')


class Test(TestCase):
    def test_load_info_file_trait(self):
        trait_info_df = load_info_file(
            logger=logger, info_file=get_path('new_ds', 'traits-lc-meta'), merge_col='Trait',
            expected_overlap_set={'Compound_287', 'Compound_287'}, reference_file='placeholder'
        )
        print(trait_info_df)

    def test_load_info_file_genes(self):
        gene_info_df = load_info_file(
            logger=logger, info_file=get_path('new_ds', 'genes-hog-info'), merge_col='Gene',
            expected_overlap_set={'N0.HOG0000000', 'N0.HOG0000001'}, reference_file='placeholder'
        )
        print(gene_info_df)

    def test_load_info_file_isolate(self):
        isolate_info_df = load_info_file(
            logger=logger, info_file=get_path('new_ds', 'isolate-meta'), merge_col='Isolate',
            expected_overlap_set={'FAM23868-i1-1.1'}, reference_file='placeholder'
        )
        print(isolate_info_df)
