from init_tests import *

from scoary.load_genes import load_genes


class Test(TestCase):
    def test_count(self):
        orig_data, binary_data = load_genes('../data/generated/Gene_presence_absence.csv', gene_data_type='gene-count')
        print(orig_data, binary_data)
        orig_data, binary_data = load_genes('../data/tetracycline/Gene_presence_absence.csv', gene_data_type='gene-count',
                                            ignore=roary_ignore)
        print(orig_data, binary_data)

    def test_list(self):
        orig_data, binary_data = load_genes(
            '../data/new_ds/Orthogroups.tsv',
            gene_data_type='gene-list:\t'
        )
        print(orig_data, binary_data)
        orig_data, binary_data = load_genes(
            '../data/new_ds/N0.tsv',
            gene_data_type='gene-list:\t',
            ignore=orthofinder_ignore
        )
        print(orig_data, binary_data)
