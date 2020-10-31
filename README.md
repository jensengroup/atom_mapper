# atom_mapper
Atom order in one molecule ("product") is made to match that in another ("reactant").

There are two main parts: `atom_mapper2D` and `atom_mapper3D`

`atom_mapper2D` is based on the idea presented in [this paper](http:/dx.doi.org/10.1021/acs.jctc.7b00764).

Note that bond orders and charges are removed when comparing fragments. I thank Leif Jacobson and  Art Bochevarov for helpful discussions.

The method breaks bonds in the reactants and products until fragments match.  The atoms in each fragments then are matched.

 `atom_mapper3D` is written by Mads Koerstz and ensures that symmetric functional groups like CH3 have the same "chirality" (when including atom labels) in reactants and products. This way the 3D coordinates from the mol object are suitable for use in reaction path interpolation methods such as NEB.
