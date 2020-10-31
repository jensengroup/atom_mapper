#
# Written by Jan Jensen 2018, 2019
#
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def change_bond_order(mol):
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)

    return mol

def change_formal_charge(mol): 
    for atom in mol.GetAtoms():    
        if atom.GetFormalCharge() != 0:
            atom.SetFormalCharge(0)
            
    return mol

def break_bonds(mol_dict,num_bonds):
    frag_smiles = []
    mol_frags = {}
    for key in mol_dict:
        mol = mol_dict[key]
        for bond in range(num_bonds):
            mol_frag = Chem.FragmentOnBonds(mol,[bond],addDummies=False)
            smiles = Chem.MolToSmiles(mol_frag,isomericSmiles=False)
            if smiles not in frag_smiles:
                frag_smiles.append(smiles)
                mol_frags[smiles] = mol_frag
                
    return frag_smiles,mol_frags
  
def get_prod_orders(matches,react_frags,prod_frags):
#
#Example: "react_frags[match]" and "prod_frags[match]" are both OC.CCN but the atom orders are
# 0-1.2-3-4 and 4-3.0-1-2.  Using GetSubstructMatch 4-3 is mapped to 0-1 and 0-1-2 is mapped to 2-3-4
#So prod_order is [4,3,0,1,2]
#
    prod_orders = []
    atom_ranks = []
    for match in matches:
        patt = react_frags[match]
        prod_order = prod_frags[match].GetSubstructMatch(patt)
        if prod_order not in prod_orders:
            prod_orders.append(prod_order)  
#            atom_ranks.append(atom_rank)
    
    return prod_orders

def atom_mapper2D(react, prod, max_bonds_cut=4):
# Remove stereochemistry in case reaction changes stereochemistry
    original_prod = Chem.Mol(prod)
    react = Chem.Mol(prod)
    Chem.rdmolops.RemoveStereochemistry(react)
    prod = Chem.Mol(prod)
    Chem.rdmolops.RemoveStereochemistry(prod)


#Remove aromaticity and work with single/double/triple bonds
    Chem.Kekulize(react,clearAromaticFlags=True)
    Chem.Kekulize(prod,clearAromaticFlags=True)
    
    react_numbonds = react.GetNumBonds()
    prod_numbonds = prod.GetNumBonds()

#Change all bonds to single bonds and remove all charges on atoms
#This allows us to compare molecule fragments based purely on connectivity
    react = change_bond_order(react)
    prod = change_bond_order(prod)
    react = change_formal_charge(react)
    prod = change_formal_charge(prod)

#recanonicalize SMILES (sanitize=False prevets hydrogen deletion)
    react_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(react),sanitize=False))
    prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(prod),sanitize=False))

#We identify active bonds based on SMILES comparisons and connect fragments to SMILES via dictionaries
    matches = []
    bonds_cut = 1
    react_frag_smiles = [react_smiles]
    prod_frag_smiles = [prod_smiles]
    react_frags = {}
    react_frags[react_smiles] = react
    prod_frags = {}
    prod_frags[prod_smiles] = prod
  
#Keep cutting bonds until the fragments match each other.
#If reactants and products have an unequal number of bonds, we cut the one with most bonds first
#until they have an equal number of bonds
    while len(matches) == 0 and bonds_cut <= max_bonds_cut:
        if react_numbonds > prod_numbonds: 
            react_frag_smiles, react_frags = break_bonds(react_frags,react_numbonds)
            react_numbonds += -1
        elif prod_numbonds > react_numbonds:
            prod_frag_smiles, prod_frags = break_bonds(prod_frags,prod_numbonds)
            prod_numbonds += -1
        else:
            react_frag_smiles, react_frags = break_bonds(react_frags,react_numbonds)
            react_numbonds += -1
            prod_frag_smiles, prod_frags = break_bonds(prod_frags,prod_numbonds)
            prod_numbonds += -1  

#Check if bond breaking produced idential set of fragments. There may be more than one way to do this
#"matches" is a list of SMILES strings describing the fragments resuling from bond cutting
        matches = list(set(react_frag_smiles) & set(prod_frag_smiles))
        bonds_cut += 1
    
    #print(matches, bonds_cut-1)

#Translate the SMILES strings to lists of atom numbers
#For example react = NCC and prod = CCN, so prod_order = [2,1,0], i.e. N is atom 0 in react and 2 in prod.
    prod_orders = get_prod_orders(matches,react_frags,prod_frags)
    #print(len(prod_orders))
        
    if len(prod_orders) == 0:
        print("no match found")
        return [i for i in range(react.GetNumAtoms())]
    if len(prod_orders) > 1:
        print("more than one match found")

    #print(list(prod_orders))

    prods_orders = reorder_prod(original_prod, prod_orders)
           
    return prods_orders

def reorder_prod(original_prod, prod_orders):
    # Reorder the atoms in the product to match that of the reactants
    prod_ordered_list = []
    for prod_order in prod_orders: 
      prod_ordered = Chem.RenumberAtoms(prod, prod_order)
      prod_ordered_list.append(prod_ordered)
      
    return prod_ordered_list

if __name__ == '__main__':
    react_smiles = "C=CC=C.C=C"
    prod_smiles = "C1CCC=CC1"

    react_smiles = "CCOOOOC"
    prod_smiles = "CCOOC.O=O"

    react_smiles = "CC=C(C)CF"
    prod_smiles = "C(C(=C)CF)C"

    #canonicalize SMILES - to get the same result independent of form of input smiles
    #react_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(react_smiles),isomericSmiles=False)
    #prod_smiles_nochiral = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles),isomericSmiles=False)
    #prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles),isomericSmiles=True)

    react = Chem.MolFromSmiles(react_smiles)
    react = Chem.AddHs(react)

    prod = Chem.MolFromSmiles(prod_smiles)
    prod = Chem.AddHs(prod)

    # Maximum number of bonds to cut when determining a match. 
    # The CPU time increases very quickly with this parameter
    max_bonds_cut = 6
    prods = atom_mapper2D(react, prod, max_bonds_cut)

    for prod in prods:
        print(Chem.MolToSmiles(prod))

