from rdkit import Chem
from rdkit.Chem import AllChem,rdMolAlign,rdMolDescriptors
import sys
from itertools import combinations, product
from collections import defaultdict


def atom_mapper(react_smiles,prod_smiles,max_bonds_cut):
# 
    prod_order = heavy_atom_mapper(react_smiles,prod_smiles_nochiral,max_bonds_cut)
    
    # Reorder the atoms in the product to match that of the reactants
    prod_ordered = Chem.RenumberAtoms(prod, prod_order)

# Overlay the most flexible molecule ("mol.sdf") on the most rigid ("template.sdf") and write sdf files
    rmsd = react_prod_rmsd(react,prod_ordered)

# Check hydrogen numbering of atoms with more than one hydrogen 
    #prod_ordered = hydrogen_atom_mapper(react,prod_ordered)
    prod_ordered = equivalent_atom_mapper(react,prod_ordered)

    return prod_ordered
    
def heavy_atom_mapper(react_smiles,prod_smiles,max_bonds_cut):
    react = Chem.MolFromSmiles(react_smiles)
    prod = Chem.MolFromSmiles(prod_smiles)

#Remove aromaticity and work with single/double/triple bonds
    Chem.Kekulize(react,clearAromaticFlags=True)
    Chem.Kekulize(prod,clearAromaticFlags=True)
        
    react = Chem.AddHs(react) 
    prod = Chem.AddHs(prod) 
    
    react_numbonds = react.GetNumBonds()
    prod_numbonds = prod.GetNumBonds()

    original_prod = prod
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
    
    print matches, bonds_cut-1

#Translate the SMILES strings to lists of atom numbers
#For example react = NCC and prod = CCN, so prod_order = [2,1,0], i.e. N is atom 0 in react and 2 in prod.
    prod_orders = get_prod_orders(matches,react_frags,prod_frags)
    print len(prod_orders)
        
    if len(prod_orders) == 0:
        print "no match found"
        return [i for i in range(react.GetNumAtoms())]
    if len(prod_orders) > 1:
        print "more than one match found"
#If more than one match is found, find the match that gives the lowest RMSD between reactant and product
        prod_order = find_best_match(react,prod,prod_orders)     
    else:
        prod_order = prod_orders[0]
    
    print list(prod_order)
           
    return list(prod_order)
        
def change_bond_order(mol):
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)

    return mol

def change_formal_charge(mol): 
    for atom in mol.GetAtoms():    
        if atom.GetFormalCharge() == 0:
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

def find_best_match(react,prod,prod_orders):
#
#If more than one match is found, find the match that gives the lowest RMSD between reactant and product     
#
    min_rmsd = 99999.9
    for prod_order in prod_orders:
        prod_ordered = Chem.RenumberAtoms(prod, prod_order)
        rmsd = react_prod_rmsd(react,prod_ordered)
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            best_prod_order = prod_order
    
    return best_prod_order

def react_prod_rmsd(react,prod, embed=True):
#Computes the RMSD between reactant and product
#The most flexible molecule ("mol") is superimposed on the most rigid ("template")

    react_rotbonds = rdMolDescriptors.CalcNumRotatableBonds(react)
    prod_rotbonds = rdMolDescriptors.CalcNumRotatableBonds(prod)
    
#If a system consists of more than one molecule then it is the most flexible
#Case where both reactant and product consists of several molecule not considered yet
#In that case arbitrarily picks product as most flexible. This migh be a problem as the 
#relative positions/orientation of fragments in template will be quite arbitrary
    if "." in Chem.MolToSmiles(react):
        mol = react
        template = prod
    elif "." in Chem.MolToSmiles(prod):
        mol = prod
        template = react
    elif react_rotbonds > prod_rotbonds:
        mol = react
        template = prod
    elif react_rotbonds < prod_rotbonds:
        mol = prod
        template = react
    else:
        mol = prod
        template = react
        
    Chem.SanitizeMol(mol)
    Chem.SanitizeMol(template)
#make a 3D model of the template and optimize w UFF
#if you got react and prod by reading in optimized geometries you want to skip these steps
    if embed:
        print "embedding"
        AllChem.EmbedMolecule(template,randomSeed=3) #remove 3 when done debugging
        AllChem.UFFOptimizeMolecule(template)

#Make a 3D model of "mol" that most closely matches that of the template an UFF optimization
#that minimizes the distances between corresponding atoms in mol and template.
    ConstrainedEmbedRP(mol, template, useTethers=False, embed=embed)
    #Chem.MolToMolBlock(mol)

    rms = float(mol.GetProp('EmbedRMS'))
    
    Chem.SDWriter("template.sdf").write(template)
    Chem.SDWriter("mol.sdf").write(mol)
    
    return rms

def ConstrainedEmbedRP(mol, core, useTethers=True, embed=True, coreConfId=-1, randomseed=2342,getForceField=AllChem.UFFGetMoleculeForceField, **kwargs):
# 
#  Copyright (C) 2006-2017  greg Landrum and Rational Discovery LLC 
# 
#   @@ All Rights Reserved @@ 
#  This file is part of the RDKit. 
#  The contents are covered by the terms of the BSD license 
#  which is included in the file license.txt, found at the root 
#  of the RDKit source tree. 
#
# Modified by Jan Jensen to work when either "mol" or "core" consists of fragments
# GetSubstructMatch only returns values for fragment with largest match and
# coordMap option doesn't work for fragments

#    match = mol.GetSubstructMatch(core)
    force_constant = 100.
    force_constant = 10.
#    match = [i for i in range(mol.GetNumAtoms())]  #jhj
    match = range(mol.GetNumAtoms()) #jhj
    if not match:
        raise ValueError("molecule doesn't match the core")
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI

    if embed:
        if "." in Chem.MolToSmiles(mol):
            ci = AllChem.EmbedMolecule(mol, randomSeed=randomseed, **kwargs)  #jhj
        else:
            ci = AllChem.EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)

        if ci < 0:
            raise ValueError('Could not embed molecule.')

    algMap = [(j, i) for i, j in enumerate(match)]

    if not useTethers:
        # clean up the conformation
        ff = getForceField(mol, confId=0)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, force_constant)
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
        ff = getForceField(mol, confId=0)
        conf = core.GetConformer()
        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0, force_constant)
        ff.Initialize()
        n = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # realign
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
    mol.SetProp('EmbedRMS', str(rms))
    return mol


def equivalent_atom_mapper(react,prod):
#inital rmsd between reactants and products
    rmsd = react_prod_rmsd(react,prod,embed=False)
    print "initial rmsd",rmsd
    
#Rank in molecular tree. Duplicates = identical atoms
#Example: CC(C)C -> [0, 3, 0, 0] meaning atom 0, 2 and 3 are identical
    atom_rank = list(Chem.CanonicalRankAtoms(prod,breakTies=False)) 
    print "atom_rank",atom_rank

# list_duplicates makes a dictionary with duplicates
# sorted(list_duplicates([0, 3, 0, 0]) -> [(0, [0, 2, 3])]
# list(itertools.permutations([0, 2, 3],2)) -> [(0, 2), (0, 3), (2, 3)]
    equivalent_atom_pairs = []
    for dummy, equivalent_atoms in list_duplicates(atom_rank):
        equivalent_atom_pairs += list(combinations(equivalent_atoms,2))
    
    print "equivalent_atom_pairs",equivalent_atom_pairs

# switch equivalent atoms pairs and save cases where rmsd goes down
    for a,b in equivalent_atom_pairs:
        order = range(prod.GetNumAtoms()) 
        order[a] = b
        order[b] = a
                
        trial_prod = Chem.RenumberAtoms(prod, order) 
        trial_rmsd = react_prod_rmsd(react,trial_prod,embed=False)
        print "rmsd, trial_rmsd",rmsd, trial_rmsd
        if trial_rmsd < rmsd:
            prod = trial_prod
            rmsd = trial_rmsd
         
    return prod

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)
            

if __name__ == '__main__':
    
    #react_smiles = "CC=C(C)CF"
    #prod_smiles = "C(C(=C)CF)C"
    
    react_smiles = "C=CC=C.C=C"
    prod_smiles = "C1CCC=CC1"
    
    #react_smiles = "C=CC=C.O=CC1=COCC1"
    #prod_smiles = "O=CC12C(CC=CC2)OCC1"
    
    #react_smiles = "C=CC=C.O=CC1=COC2=CC=CC=C21"
    #prod_smiles = "O=CC12C(CC=CC2)OC3=CC=CC=C31"
    
    #react_smiles = "C/C(C)=C/CC/C(C)=C/CC/C(C)=C/C"
    #prod_smiles = "CC1(C)CCC[C@@]2(C)[C@@]1([H])CCC([C@@H]2C)=C"
    #prod_smiles = "CC1(C)CCCC2(C)C1CCC(C2C)=C"
    
    #react_smiles = "N#CC1(C#N)C(C2=CC=CC=C2)=CC3=CC=CC=CC31"
    #prod_smiles = "N#C/C(C#N)=C(C1=CC=CC=C1)\C=C2C=CC=CC=C2"
    
    #react_smiles = "N#CC1(C#N)C=CC3=CC=CC=CC31"
    #prod_smiles = "N#C/C(C#N)=C\C=C2C=CC=CC=C2"

#canonicalize SMILES - to get the same result independent of form of input smiles
    react_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(react_smiles),isomericSmiles=False)
    prod_smiles_nochiral = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles),isomericSmiles=False)
    prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles),isomericSmiles=True)
    
    react = Chem.MolFromSmiles(react_smiles)
    react = Chem.AddHs(react)

    prod = Chem.MolFromSmiles(prod_smiles)
    prod = Chem.AddHs(prod)

# Maximum number of bonds to cut when determining a match. 
# The CPU time increases very quickly with this parameter
    max_bonds_cut = 6
    prod_ordered = atom_mapper(react_smiles,prod_smiles_nochiral,max_bonds_cut)
  
    rmsd = react_prod_rmsd(react,prod_ordered,embed=False)
    print "final rmsd", rmsd
