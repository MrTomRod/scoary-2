from re import compile

BRANCH_LENTGHS = compile(r':[0-9]+(.[0-9]+)?')


class NewickParserException(Exception):
    pass


def parse_newick(newick_string: str) -> []:
    """
    A simple function to parse Newick strings to list tree format.

    Example:
        >>> parse_newick('(A,(B,C))D;')
        [['A', ['B', 'C']], 'D']

    Limitations:
        - Only binary trees are supported
        - Distances are ignored
        - All labels must be named
        - NHX format not supported

    :param newick_string: Phylogenetic tree in newick format
    :return: Phylogenetic tree in list format
    """

    # strip and remove branch lengths
    newick_string = newick_string.strip()
    newick_string = BRANCH_LENTGHS.sub(string=newick_string, repl='')

    # sanity check
    if not newick_string.endswith(';'):
        raise NewickParserException(f'Newick string does not end in semicolon! {newick_string=}')

    def find_corresponding_closing(string: str) -> int:
        n_opening = 0
        n_closing = 0
        for i, char in enumerate(string):
            if char == '(':
                n_opening += 1
                continue
            if char == ')':
                n_closing += 1
                if n_closing == n_opening:
                    return i

        raise NewickParserException(f'Could not find corresponding closing bracket in {string=}! {newick_string=}')

    def split_node(string: str) -> (str, str):
        if ',' in string:
            first_comma = string.index(',')
            if '(' in string:
                first_bracket = string.index('(')
                if first_bracket < first_comma:
                    return string[0:first_bracket], string[first_bracket:]

            return string[0:first_comma], string[first_comma + 1:]
        else:
            if '(' not in string:
                raise NewickParserException(f'Could not find separators "," or "(" in {string=}! {newick_string=}')
            first_bracket = string.index('(')
            return string[0:first_bracket], string[first_bracket:]

    def parse_leaf(string: str) -> str:
        string = string.strip('"\' ')
        if len(string) == 0:
            raise NewickParserException(f'Leaf with no label: {string=}! {newick_string=}')
        return string

    def parse_recursive(string: str) -> list | str:
        string = string.strip()

        # leaf
        if not ',' in string and not '(' in string and not ')' in string:
            return parse_leaf(string)

        # remove enclosing brackets
        if string.startswith('(') and string.endswith(')') \
                and find_corresponding_closing(string) == len(string) - 1:
            string = string[1:-1]

        # parse node
        if string.startswith('('):
            closing_idx = find_corresponding_closing(string)
            left = string[1:closing_idx]
            right = string[closing_idx + 1:].lstrip(',')
        else:
            left, right = split_node(string)

        left = parse_recursive(left)
        right = parse_recursive(right)

        return [left, right]

    return parse_recursive(newick_string[:-1])  # remove semicolon
