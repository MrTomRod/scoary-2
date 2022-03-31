from unittest import TestCase

import pandas as pd

from init_tests import *

from scoary.load_genes import load_genes


class Test(TestCase):
    def test_count(self):
        orig_data, binary_data = load_genes(get_path('generated', 'genes'), gene_data_type='gene-count')
        print(orig_data, binary_data)
        orig_data, binary_data = load_genes(get_path('tetracycline', 'genes'), gene_data_type='gene-count',
                                            ignore=tetr_ignore)
        print(orig_data, binary_data)

    def test_list(self):
        orig_data, binary_data = load_genes(get_path('new_ds', 'genes-og'), gene_data_type='gene-list')
        print(orig_data, binary_data)
        orig_data, binary_data = load_genes(get_path('new_ds', 'genes-hog'), gene_data_type='gene-list',
                                            ignore=tetr_ignore)
        print(orig_data, binary_data)

    def test_list_annotations(self):
        orig_data, binary_data, info = load_genes(
            get_path('new_ds', 'genes-hog'),
            gene_data_type='gene-list',
            gene_info=get_path('new_ds', 'genes-hog-info')
        )
        print(orig_data, binary_data, info)
