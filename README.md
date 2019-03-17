# atom_mapper
Atom order in one molecule ("product") is made to match that in another ("reactant").

The code suggests a reaction coordinate for a constrained scan.

Based on the idea presented in [this paper](http:/dx.doi.org/10.1021/acs.jctc.7b00764).

Note that bond orders and charges are removed when comparing fragments. I thank Leif Jacobson and  Art Bochevarov for helpful discussions.

The method breaks bonds in the reactants and products until fragments match.  The atoms in each fragments then are matched.

Note that if you generate 3D coordinates from the mol object, then symmetric functional groups like CH3 may have different "chirality" (when including atom labels) in reactants and products. So these coordinates may not be suitable for interpolation. 
