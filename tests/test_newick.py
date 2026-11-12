import pytest

from scoary.newick import parse_newick, NewickParserException
from conftest import TETRACYCLINE


@pytest.mark.parametrize('newick, expected', [
    ('(A,B);', ['A', 'B']),
    ('(A,(B,C))D;', [['A', ['B', 'C']], 'D']),
    ('((A,B),(C,D));', [['A', 'B'], ['C', 'D']]),
    ('(A,(B,(D,G)));', ['A', ['B', ['D', 'G']]]),
    # branch lengths, incl. scientific notation, are dropped
    ('(A:0.1,(B:0.2,C:3e-05):0.4);', ['A', ['B', 'C']]),
    # support values after closing brackets are dropped
    ('(FAM18356-i1-1.1:0.289313,(FAM19471-i1-1.1:0.210374,(FAM23169-i1-1.1:0.0764835,'
     'FAM1079-i1-1.1:0.0778624)0.977444:0.115601));',
     ['FAM18356-i1-1.1', ['FAM19471-i1-1.1', ['FAM23169-i1-1.1', 'FAM1079-i1-1.1']]]),
    # labels with spaces, dashes, underscores, leading digits
    ('(hodor one,(hodor-two,(hodor_3,0hodor 4)));',
     ['hodor one', ['hodor-two', ['hodor_3', '0hodor 4']]]),
    # quoted labels, zero-length branches, support on root
    ("((('A':1.48e-06,('B':0,'C':0)0.98:0)1:0,'D':0)0.94:2.32e-06,'E':0)1;",
     [[['A', ['B', 'C']], 'D'], 'E']),
    # surrounding whitespace / trailing newline
    ('  (A,B);\n', ['A', 'B']),
])
def test_parse_newick(newick, expected):
    assert parse_newick(newick) == expected


@pytest.mark.parametrize('bad', [
    '(A,B)',  # missing semicolon
    '(A,);',  # empty leaf
    '((A,B);',  # unbalanced brackets
])
def test_parse_newick_rejects_malformed(bad):
    with pytest.raises(NewickParserException):
        parse_newick(bad)


def test_parse_example_tree():
    tree = parse_newick((TETRACYCLINE / 'ExampleTree.nwk').read_text())

    def leaves(t):
        return [t] if isinstance(t, str) else leaves(t[0]) + leaves(t[1])

    labels = leaves(tree)
    assert len(labels) == len(set(labels)) == 100
    assert all(l.startswith('Isolate_') for l in labels)
