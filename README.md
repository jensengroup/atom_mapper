# atom_mapper
Atom order in one molecule is made to match that in another

Based on the idea presented in this paper: 10.1021/acs.jctc.7b00764

The method breaks bonds in the reactants and products until fragments match.  The the atoms in each fragments then are matched.

Example: fragmenting CC=C(C)CF and C(C(=C)CF)C leads to the identical fragment C.C=C(C)CF after breaking 1 bond each, where the atoms in "C" and "C=C(C)CF" can be easily matched because they are identical.

TODO:
1. The fragmentation can lead to more than one of the same fragment type. E.g. for the prototypical Diels-Alder reaction C=CC=C.C=C>>C1CCCC=C1 two fragmentations pattern match: C.C.C=C.CC' and 'C.C.C.C=CC'.  Also there is more than one "C" fragment.
The paper suggests trying all permutations to see which one results in the lowest RMSD

2. Determining atom order before adding Hs, which is more efficient, seems to work, but hasn't been tested very much.
