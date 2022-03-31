class PhyloTree:
    """
    A class that represents a binary tree. Phylotrees can be nested.
    They can also contain tips at left, right or both nodes. Stores the
    max number of paths under each of the following 5 conditions:
    1. Free path to AB.
    2. Free path to Ab.
    3. Free path to aB.
    4. Free path to ab.
    5. No free path.
    """

    def __init__(
            self,
            leftnode,  # ['strainA', 'strainB'] or 'strainA'
            rightnode,  # ['strainA', 'strainB'] or 'strainA'
            GTC  # {'strain1': 'ab', 'strain2': 'Ab', 'strain3': 'ab'}
    ):
        """
        Constructs a phylotree and links it to its left and right nodes
        """
        # First check if left and right are tips (type = str). Then we
        # need to look up the gene-trait combination value at that node.
        # Also check if left and right nodes are PhyloTrees. If not,
        # recursively create PhyloTrees.

        if type(leftnode) is str:
            self.leftnode = Tip(GTC[leftnode], leftnode)
        elif isinstance(leftnode, PhyloTree):
            self.leftnode = leftnode
        else:
            self.leftnode = PhyloTree(leftnode=leftnode[0],
                                      rightnode=leftnode[1],
                                      GTC=GTC)

        if type(rightnode) is str:
            self.rightnode = Tip(GTC[rightnode], rightnode)
        elif isinstance(rightnode, PhyloTree):
            self.rightnode = rightnode
        else:
            self.rightnode = PhyloTree(leftnode=rightnode[0],
                                       rightnode=rightnode[1],
                                       GTC=GTC)

        # Initialize the max number of paths. Set to -1 meaning they
        # cannot be reached
        self.maxvalues = \
            {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.max_propairs = \
            {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.max_antipairs = \
            {"AB": -1, "Ab": -1, "aB": -1, "ab": -1, "0": -1}
        self.calculate_max()
        self.max_contrasting_pairs = max(self.maxvalues.values())
        self.max_contrasting_propairs = max(self.max_propairs.values())
        self.max_contrasting_antipairs = max(self.max_antipairs.values())

    def __str__(self):
        return f'({self.leftnode}, {self.rightnode})'

    def _to_newick(self):
        return f'({self.leftnode._to_newick()}, {self.rightnode._to_newick()})'

    def to_newick(self):
        return f'({self.leftnode._to_newick()}, {self.rightnode._to_newick()});'

    def to_label_to_gene(self):
        return self.leftnode.to_label_to_gene() | self.rightnode.to_label_to_gene()

    def to_label_to_trait(self):
        return self.leftnode.to_label_to_trait() | self.rightnode.to_label_to_trait()

    def calculate_max(self):
        """
        A method for calculating the max number of pairings under the 5
        conditions
        """

        for condition in ["AB", "Ab", "aB", "ab"]:
            pairings = self.calculate_max_condition(condition)  # {"Total": int, "Pro": int, "Anti": int}
            self.maxvalues[condition] = pairings["Total"]
            self.max_propairs[condition] = pairings["Pro"]
            self.max_antipairs[condition] = pairings["Anti"]
        pairings = self.calculate_max_nofree()
        self.maxvalues["0"] = pairings["Total"]
        self.max_propairs["0"] = pairings["Pro"]
        self.max_antipairs["0"] = pairings["Anti"]

    def calculate_max_condition(self, condition):
        """
        When passed for example 'AB', finds out the 9 distinct
        conditions and calculates the max. Here, the possibilities are
        illustrated when passed 'AB' (eg. a free path to AB exists):

        Left        Right
        AB          No free
        AB          Ab
        AB          aB
        AB          ab
        AB          AB
        No free     AB
        Ab          AB
        aB          AB
        ab          AB

        """
        Possible_conditions = {"AB", "Ab", "aB", "ab"}
        Possible_conditions.remove(condition)
        # Now we have a list of the elements that are NOT condition:
        other_conditions = list(Possible_conditions)
        max_pairs_1 = -1
        max_pairs_2 = -1
        max_pairs_3 = -1
        max_pairs_4 = -1
        max_pairs_5 = -1
        max_pairs_6 = -1
        max_pairs_7 = -1
        max_pairs_8 = -1
        max_pairs_9 = -1

        max_propairs_1 = -1
        max_propairs_2 = -1
        max_propairs_3 = -1
        max_propairs_4 = -1
        max_propairs_5 = -1
        max_propairs_6 = -1
        max_propairs_7 = -1
        max_propairs_8 = -1
        max_propairs_9 = -1

        max_antipairs_1 = -1
        max_antipairs_2 = -1
        max_antipairs_3 = -1
        max_antipairs_4 = -1
        max_antipairs_5 = -1
        max_antipairs_6 = -1
        max_antipairs_7 = -1
        max_antipairs_8 = -1
        max_antipairs_9 = -1

        l_maxvals = self.leftnode.maxvalues
        r_maxvals = self.rightnode.maxvalues

        l_maxpros = self.leftnode.max_propairs
        r_maxpros = self.rightnode.max_propairs

        l_maxantis = self.leftnode.max_antipairs
        r_maxantis = self.rightnode.max_antipairs

        # special case: comparison with nonfree
        if l_maxvals[condition] > -1 and r_maxvals["0"] > -1:
            max_pairs_1 = l_maxvals[condition] + r_maxvals["0"]
            max_propairs_1 = l_maxpros[condition] + r_maxpros["0"]
            max_antipairs_1 = l_maxantis[condition] + r_maxantis["0"]

        # comparison with all others
        if l_maxvals[condition] > -1 and r_maxvals[other_conditions[0]] > -1:
            max_pairs_2 = l_maxvals[condition] + r_maxvals[other_conditions[0]]
            max_propairs_2 = l_maxpros[condition] + r_maxpros[other_conditions[0]]
            max_antipairs_2 = l_maxantis[condition] + r_maxantis[other_conditions[0]]

        if l_maxvals[condition] > -1 and r_maxvals[other_conditions[1]] > -1:
            max_pairs_3 = l_maxvals[condition] + r_maxvals[other_conditions[1]]
            max_propairs_3 = l_maxpros[condition] + r_maxpros[other_conditions[1]]
            max_antipairs_3 = l_maxantis[condition] + r_maxantis[other_conditions[1]]

        if l_maxvals[condition] > -1 and r_maxvals[other_conditions[2]] > -1:
            max_pairs_4 = l_maxvals[condition] + r_maxvals[other_conditions[2]]
            max_propairs_4 = l_maxpros[condition] + r_maxpros[other_conditions[2]]
            max_antipairs_4 = l_maxantis[condition] + r_maxantis[other_conditions[2]]

        # special case: comparison with self
        if l_maxvals[condition] > -1 and r_maxvals[condition] > -1:
            max_pairs_5 = l_maxvals[condition] + r_maxvals[condition]
            max_propairs_5 = l_maxpros[condition] + r_maxpros[condition]
            max_antipairs_5 = l_maxantis[condition] + r_maxantis[condition]

        # comparison with all others
        if l_maxvals["0"] > -1 and r_maxvals[condition] > -1:
            max_pairs_6 = l_maxvals["0"] + r_maxvals[condition]
            max_propairs_6 = l_maxpros["0"] + r_maxpros[condition]
            max_antipairs_6 = l_maxantis["0"] + r_maxantis[condition]

        if l_maxvals[other_conditions[0]] > -1 and r_maxvals[condition] > -1:
            max_pairs_7 = l_maxvals[other_conditions[0]] + r_maxvals[condition]
            max_propairs_7 = l_maxpros[other_conditions[0]] + r_maxpros[condition]
            max_antipairs_7 = l_maxantis[other_conditions[0]] + r_maxantis[condition]

        if l_maxvals[other_conditions[1]] > -1 and r_maxvals[condition] > -1:
            max_pairs_8 = l_maxvals[other_conditions[1]] + r_maxvals[condition]
            max_propairs_8 = l_maxpros[other_conditions[1]] + r_maxpros[condition]
            max_antipairs_8 = l_maxantis[other_conditions[1]] + r_maxantis[condition]

        if l_maxvals[other_conditions[2]] > -1 and r_maxvals[condition] > -1:
            max_pairs_9 = l_maxvals[other_conditions[2]] + r_maxvals[condition]
            max_propairs_9 = l_maxpros[other_conditions[2]] + r_maxpros[condition]
            max_antipairs_9 = l_maxantis[other_conditions[2]] + r_maxantis[condition]

        max_pairs = max(max_pairs_1, max_pairs_2, max_pairs_3, max_pairs_4, max_pairs_5,
                        max_pairs_6, max_pairs_7, max_pairs_8, max_pairs_9)

        # Calculate maximum number of propairs, given a maxmimum number
        # of pairs
        max_propairs = -1
        if max_pairs == max_pairs_1:
            max_propairs = max(max_propairs, max_propairs_1)
        if max_pairs == max_pairs_2:
            max_propairs = max(max_propairs, max_propairs_2)
        if max_pairs == max_pairs_3:
            max_propairs = max(max_propairs, max_propairs_3)
        if max_pairs == max_pairs_4:
            max_propairs = max(max_propairs, max_propairs_4)
        if max_pairs == max_pairs_5:
            max_propairs = max(max_propairs, max_propairs_5)
        if max_pairs == max_pairs_6:
            max_propairs = max(max_propairs, max_propairs_6)
        if max_pairs == max_pairs_7:
            max_propairs = max(max_propairs, max_propairs_7)
        if max_pairs == max_pairs_8:
            max_propairs = max(max_propairs, max_propairs_8)
        if max_pairs == max_pairs_9:
            max_propairs = max(max_propairs, max_propairs_9)

        # Calculate maximum number of antipairs, given a maxmimum number
        # of pairs
        max_antipairs = -1
        if max_pairs == max_pairs_1:
            max_antipairs = max(max_antipairs, max_antipairs_1)
        if max_pairs == max_pairs_2:
            max_antipairs = max(max_antipairs, max_antipairs_2)
        if max_pairs == max_pairs_3:
            max_antipairs = max(max_antipairs, max_antipairs_3)
        if max_pairs == max_pairs_4:
            max_antipairs = max(max_antipairs, max_antipairs_4)
        if max_pairs == max_pairs_5:
            max_antipairs = max(max_antipairs, max_antipairs_5)
        if max_pairs == max_pairs_6:
            max_antipairs = max(max_antipairs, max_antipairs_6)
        if max_pairs == max_pairs_7:
            max_antipairs = max(max_antipairs, max_antipairs_7)
        if max_pairs == max_pairs_8:
            max_antipairs = max(max_antipairs, max_antipairs_8)
        if max_pairs == max_pairs_9:
            max_antipairs = max(max_antipairs, max_antipairs_9)

        return {"Total": max_pairs, "Pro": max_propairs, "Anti": max_antipairs}

    def calculate_max_nofree(self):
        """
        Under the condition of no free paths, only 5 distinct
        possibilities exits. (No free paths requires either that there
        are no free paths in the left or right nodes, or that a new pair
        is formed across the root from the left and right nodes)

        Left        Right       Result
        No free     No free
        AB          ab          +1 pair (pro)
        ab          AB          +1 pair (pro)
        Ab          aB          +1 pair (anti)
        aB          Ab          +1 pair (anti)
        """
        # No free pairs in either:
        max_pairs_nofree = -1
        max_pairs_1100 = -1
        max_pairs_0011 = -1
        max_pairs_1001 = -1
        max_pairs_0110 = -1

        max_propairs_nofree = -1
        max_propairs_1100 = -1
        max_propairs_0011 = -1
        max_propairs_1001 = -1
        max_propairs_0110 = -1

        max_antipairs_nofree = -1
        max_antipairs_1100 = -1
        max_antipairs_0011 = -1
        max_antipairs_1001 = -1
        max_antipairs_0110 = -1

        l_maxvals = self.leftnode.maxvalues
        r_maxvals = self.rightnode.maxvalues

        l_maxpros = self.leftnode.max_propairs
        r_maxpros = self.rightnode.max_propairs

        l_maxantis = self.leftnode.max_antipairs
        r_maxantis = self.rightnode.max_antipairs

        if l_maxvals["0"] > -1 and r_maxvals["0"] > -1:
            max_pairs_nofree = l_maxvals["0"] + r_maxvals["0"]
            max_propairs_nofree = l_maxpros["0"] + r_maxpros["0"]
            max_antipairs_nofree = l_maxantis["0"] + r_maxantis["0"]

        if l_maxvals["AB"] > -1 and r_maxvals["ab"] > -1:
            max_pairs_1100 = l_maxvals["AB"] + r_maxvals["ab"] + 1
            max_propairs_1100 = l_maxpros["AB"] + r_maxpros["ab"] + 1
            max_antipairs_1100 = l_maxantis["AB"] + r_maxantis["ab"]

        if l_maxvals["ab"] > -1 and r_maxvals["AB"] > -1:
            max_pairs_0011 = l_maxvals["ab"] + r_maxvals["AB"] + 1
            max_propairs_0011 = l_maxpros["ab"] + r_maxpros["AB"] + 1
            max_antipairs_0011 = l_maxantis["ab"] + r_maxantis["AB"]

        if l_maxvals["Ab"] > -1 and r_maxvals["aB"] > -1:
            max_pairs_1001 = l_maxvals["Ab"] + r_maxvals["aB"] + 1
            max_propairs_1001 = l_maxpros["Ab"] + r_maxpros["aB"]
            max_antipairs_1001 = l_maxantis["Ab"] + r_maxantis["aB"] + 1

        if l_maxvals["aB"] > -1 and r_maxvals["Ab"] > -1:
            max_pairs_0110 = l_maxvals["aB"] + r_maxvals["Ab"] + 1
            max_propairs_0110 = l_maxpros["aB"] + r_maxpros["Ab"]
            max_antipairs_0110 = l_maxantis["aB"] + r_maxantis["Ab"] + 1

        max_pairs = max(max_pairs_nofree, max_pairs_1100,
                        max_pairs_0011, max_pairs_1001, max_pairs_0110)

        # Calculate max number of propairs
        max_propairs = -1  # Max_propairs can never go below -1
        if max_pairs == max_pairs_nofree:
            max_propairs = max(max_propairs, max_propairs_nofree)
        if max_pairs == max_pairs_1100:
            max_propairs = max(max_propairs, max_propairs_1100)
        if max_pairs == max_pairs_0011:
            max_propairs = max(max_propairs, max_propairs_0011)
        if max_pairs == max_pairs_1001:
            max_propairs = max(max_propairs, max_propairs_1001)
        if max_pairs == max_pairs_0110:
            max_propairs = max(max_propairs, max_propairs_0110)

        # Calculate max number of antipairs
        max_antipairs = -1  # Max_antipairs can never go below -1
        if max_pairs == max_pairs_nofree:
            max_antipairs = max(max_antipairs, max_antipairs_nofree)
        if max_pairs == max_pairs_1100:
            max_antipairs = max(max_antipairs, max_antipairs_1100)
        if max_pairs == max_pairs_0011:
            max_antipairs = max(max_antipairs, max_antipairs_0011)
        if max_pairs == max_pairs_1001:
            max_antipairs = max(max_antipairs, max_antipairs_1001)
        if max_pairs == max_pairs_0110:
            max_antipairs = max(max_antipairs, max_antipairs_0110)

        return {"Total": max_pairs, "Pro": max_propairs, "Anti": max_antipairs}


class Tip:
    """
    A class that references a single tip, which can only be AB, Ab, aB,
    ab, A- or a-
    """

    def __init__(self, tipvalue, tipname: str):
        """
        Sets up the tip
        """
        self.tipvalue = tipvalue  # AB
        self.tipname = tipname  # strainA
        self.maxvalues = {}
        for condition in ["AB", "Ab", "aB", "ab", "0"]:
            if condition == tipvalue:
                self.maxvalues[condition] = 0
            else:
                self.maxvalues[condition] = -1
        self.max_propairs = {k: v for (k, v) in self.maxvalues.items()}
        self.max_antipairs = {k: v for (k, v) in self.maxvalues.items()}

    def __str__(self):
        tv1, tv2 = self.tipvalue
        assert tv1 in 'aA' and tv2 in 'bB'
        tv1 = '1' if tv1 == 'A' else '0'
        tv2 = '1' if tv2 == 'B' else '0'
        return f'{tv1}{tv2}_{self.tipname}'

    def _to_newick(self):
        return self.tipname

    def to_label_to_gene(self):
        tv1, tv2 = self.tipvalue
        assert tv1 in 'aA' and tv2 in 'bB'
        return {self.tipname: tv1 == 'A'}

    def to_label_to_trait(self):
        tv1, tv2 = self.tipvalue
        assert tv1 in 'aA' and tv2 in 'bB'
        return {self.tipname: tv2 == 'B'}


def convert_upgma_to_phylotree(tree, GTC):
    """
    A method that converts the upgma tree (in nested list form) to a
    PhyloTree. It also needs the status (AB, Ab, aB or ab) of all
    strains that are to be included in the analysis. As of 1.6.0 this
    function is also called from permutations
    """

    # TRAVERSING TREE: For each binary division - go to left until hit
    # tip. Then go back
    MyPhyloTree = PhyloTree(leftnode=tree[0],
                            rightnode=tree[1],
                            GTC=GTC)

    return MyPhyloTree, {"Total": MyPhyloTree.max_contrasting_pairs,
                         "Pro": MyPhyloTree.max_contrasting_propairs,
                         "Anti": MyPhyloTree.max_contrasting_antipairs}


import random


def permute_gtc(GTC):
    """
    Returns a permutation of the gene-trait combination dic where trait
    labels have been swapped

    ISOLATE KEEP GENES, GET PERMUTED TRAITS. TOTAL NUMBER OF TRAITS REMAINS
    """
    # Set of trait labels to distribute
    trait_labels = [s[-1] for s in GTC.values()]
    # Shuffle list of labels
    random.shuffle(trait_labels)
    for isolate in GTC:
        # Assign the latter character of GTC to isolates but keep the
        # gene status intact
        GTC[isolate] = str(GTC[isolate][0] + trait_labels.pop())
    return GTC
