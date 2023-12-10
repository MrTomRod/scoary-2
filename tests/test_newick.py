from init_tests import *

from scoary.newick import parse_newick


class Test(TestCase):
    def test_newick(self):
        test_data = [
            ('(A,(C,D));', ['A', ['C', 'D']]),
            ('(A, (C,D));', ['A', ['C', 'D']]),
            ('(A(C,D));', ['A', ['C', 'D']]),
            ('(A(C, D));', ['A', ['C', 'D']]),
            ('A,(C,D);', ['A', ['C', 'D']]),
            ('((A,B),(C,D));', [['A', 'B'], ['C', 'D']]),
            ('(A,B),(C,D);', [['A', 'B'], ['C', 'D']]),
            ('(A,B),(C,D);', [['A', 'B'], ['C', 'D']]),
            ('(A,B)(C,D);', [['A', 'B'], ['C', 'D']]),
            ('(C,D)E;', [['C', 'D'], 'E']),
            ('(C,D),E;', [['C', 'D'], 'E']),
            ('(A,(C,D))F;', [['A', ['C', 'D']], 'F']),
            ('(  A  , (  C  ,  D ) ) F ;', [['A', ['C', 'D']], 'F']),
            ('(A(C,D))F;', [['A', ['C', 'D']], 'F']),
            ('(A:0.1,(C:0.3,D:0.4):0.5);', ['A', ['C', 'D']]),
            ('(A:0.1,(C:0.3,D:0.4))F;', [['A', ['C', 'D']], 'F']),
            ('((B:0.2,(C:0.3,D:0.4))F:0.1)A;', [[['B', ['C', 'D']], 'F'], 'A']),
            ('(A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);', ['A', [['B', [['D', 'G'], 'E']], 'C']]),
            ('(A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);', ['A', [['B', [['D', 'G'], 'E']], 'C']]),
            ('(A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)));', ['A', ['B', ['D', 'G']]]),
            ('(A:0.350596,(B:0.728431,(D:0.609498,G:0.125729):0.642905):0.567737);', ['A', ['B', ['D', 'G']]]),
            ('(A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E)C);', ['A', [['B', [['D', 'G'], 'E']], 'C']]),
            ('(A,(B,(D,G)E)C);', ['A', [['B', [['D', 'G'], 'E']], 'C']]),
            ('(A,(B,(D,G)));', ['A', ['B', ['D', 'G']]]),
            ('(hodor one,(hodor-two,(hodor_3,0hodor 4)));', ['hodor one', ['hodor-two', ['hodor_3', '0hodor 4']]]),
        ]

        for n, expected_result in test_data:
            res = parse_newick(n)
            self.assertEqual(res, expected_result)
            print()

    def test_tetracycline(self):
        with open('../data/tetracycline/ExampleTree.nwk') as f:
            newick = f.read()
        expected_result = get_json('../data/tetracycline/expected_result.json')['as_list']
        res = parse_newick(newick)
        self.assertEqual(res, expected_result)
