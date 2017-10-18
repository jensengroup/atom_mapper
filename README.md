# atom_mapper
Atom order in one molecule is made to match that in another

Based on the idea presented in this paper: 10.1021/acs.jctc.7b00764
Note that bond orders and charges are removed which is not mentioned in the paper.

The method breaks bonds in the reactants and products until fragments match.  The atoms in each fragments then are matched.

If there is more than one match, the match that produces the lowest RMSD between reactants and products is chosen.

The set of reactant and producut coordinates with the lowest RMSD are written out as "mol.sdf" and "template.sdf" where template.sdf is the molecule judged to be most rigid.

In the case where both reactants and products consists of more than one molecule, the product is arbitrarily picked as most flexible. This migh be a problem as the relative positions/orientation of fragments in template will be quite arbitrary
