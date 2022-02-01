import pandas as pd

FISHER_COMBINATIONS = ['abcd', 'acbd', 'badc', 'bdac', 'cadb', 'cdab', 'dbca', 'dcba']


def fisher_id(a: int, b: int, c: int, d: int) -> (int, int, int, int):
    """
    Eight contingency tables always give the same pvalue: ['abcd', 'acbd', 'badc', 'bdac', 'cadb', 'cdab', 'dbca', 'dcba']

    Compute and save only one version.
    """
    vals = {'a': a, 'b': b, 'c': c, 'd': d}
    equivalent_combinations = [tuple(vals[letter] for letter in combination) for combination in FISHER_COMBINATIONS]
    return sorted(equivalent_combinations, key=lambda comb: (comb[0], comb[1], comb[2], comb[3]))[0]
