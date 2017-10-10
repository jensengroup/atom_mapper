from rdkit import Chem
from rdkit.Chem import AllChem

def atom_mapper(react_smiles,prod_smiles,max_bonds_cut):
    react = Chem.MolFromSmiles(react_smiles)
    prod = Chem.MolFromSmiles(prod_smiles)
    
    match = []
    bonds_cut = 1
    react_frag_smiles = [react_smiles]
    prod_frag_smiles = [prod_smiles]
    react_frags = {}
    react_frags[react_smiles] = react
    prod_frags = {}
    prod_frags[react_smiles] = prod
    react_numbonds = react.GetNumBonds()
    prod_numbonds = prod.GetNumBonds()
    
    while len(match) == 0 and bonds_cut <= max_bonds_cut:
        if react_numbonds >= prod_numbonds: 
            react_frag_smiles, react_frags = break_bonds(react_frags,react_numbonds)
            react_numbonds += -1
        if prod_numbonds >= react_numbonds:
            prod_frag_smiles, prod_frags = break_bonds(prod_frags,prod_numbonds)
            prod_numbonds += -1
                
        match = list(set(react_frag_smiles) & set(prod_frag_smiles))
        bonds_cut += 1
    
    print match
        
    if len(match) == 0:
        print "no match found"
        return prod_smiles
    if len(match) > 1:
        print "more than one match found"
        return prod_smiles        

    patt = react_frags[match[0]]
    prod_order = prod_frags[match[0]].GetSubstructMatches(patt)
    
    print prod_order[0]
    
    Chem.SanitizeMol(prod)
    prod_ordered = Chem.RenumberAtoms(prod, list(prod_order[0]))
    prod_ordered_smiles = Chem.MolToSmiles(prod_ordered,canonical=False)

    return prod_ordered_smiles
        
def break_bonds(mol_dict,num_bonds):
    frag_smiles = []
    mol_frags = {}
    for key in mol_dict:
        mol = mol_dict[key]
        for bond in range(num_bonds):
            mol_frag = Chem.FragmentOnBonds(mol,[bond],addDummies=False)
            smiles = Chem.MolToSmiles(mol_frag,True)
            if smiles not in frag_smiles:
                frag_smiles.append(smiles)
                mol_frags[smiles] = mol_frag
                
    return frag_smiles,mol_frags


if __name__ == '__main__':
    max_bonds_cut = 6
    react_smiles = "CC=C(C)CF"
    prod_smiles = "C(C(=C)CF)C"
    
    react_smiles = "C=CC=C.C=C"
    prod_smiles = "C1CCC=CC1"
    
    react_smiles = "C/C(C)=C/CC/C(C)=C/CC/C(C)=C/C"
    prod_smiles = "CC1(C)CCC[C@@]2(C)[C@@]1([H])CCC([C@@H]2C)=C"
    #prod_smiles = "CC1(C)CCCC2(C)C1CCC(C2C)=C"

#canonicalize SMILES
    react_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(react_smiles))
    prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles))

    prod_ordered_smiles = atom_mapper(react_smiles,prod_smiles,max_bonds_cut)

    print prod_ordered_smiles
