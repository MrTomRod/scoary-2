from .init_tests import *

from scoary import ScoaryTree, pick_single, print_tree, pick

from scoary.permutations import create_permuted_df, permute_picking


def test_permutations(tree: list, label_to_trait_a, label_to_trait_b, n_permut):
    max_contr, max_suppo, max_oppos, best, worst = pick(
        tree=tree,
        label_to_trait_a=label_to_trait_a,
        trait_b_df=pd.DataFrame(label_to_trait_b, index=['gene']),
        calc_pvals=True
    )
    print_tree(
        ScoaryTree.from_list(tree),
        label_to_trait_a, label_to_trait_b
    )
    is_positively_correlated = max_suppo >= max_oppos
    n_pos, n_neg = sum(label_to_trait_b.values()), len(label_to_trait_b)
    n_positive = n_pos if is_positively_correlated else n_neg

    estimator = (max_suppo if is_positively_correlated else max_oppos) / max_contr
    print(f'{max_contr=}\n{max_suppo=}\n{max_oppos=}\n{best=}\n{worst=}\n{estimator=}')

    print('Calculating permutatons... p-value=', end='')
    permuted_df = create_permuted_df(
        labels=[f'i{i}' for i in range(1, 17)], n_positive=n_positive,
        n_permut=n_permut, random_state=42
    )
    max_contr, max_suppo, max_oppos = pick(
        tree=tree, label_to_trait_a=label_to_trait_a,
        trait_b_df=permuted_df, calc_pvals=False
    )

    permuted_estimators = max_suppo / max_contr

    pval = ((permuted_estimators >= estimator).sum() + 1) / (n_permut + 1)

    print(pval)


class Test(TestCase):
    def test_bad(self, n_permut=3000):
        tree = [[[['i1', 'i2'], ['i3', 'i4']], [['i5', 'i6'], ['i7', 'i8']]],
                [[['i9', 'i10'], ['i11', 'i12']], [['i13', 'i14'], ['i15', 'i16']]]]
        label_to_trait_a = {f'i{i}': bool(i % 2) for i in range(1, 17)}
        label_to_trait_b = label_to_trait_a.copy()
        test_permutations(tree, label_to_trait_a, label_to_trait_b, n_permut)

    def test_good(self, n_permut=3000):
        tree = [[[['i1', 'i2'], ['i3', 'i4']], [['i5', 'i6'], ['i7', 'i8']]],
                [[['i9', 'i10'], ['i11', 'i12']], [['i13', 'i14'], ['i15', 'i16']]]]
        label_to_trait_a = {f'i{i}': bool(i < 9) for i in range(1, 17)}
        label_to_trait_b = label_to_trait_a.copy()
        test_permutations(tree, label_to_trait_a, label_to_trait_b, n_permut)
