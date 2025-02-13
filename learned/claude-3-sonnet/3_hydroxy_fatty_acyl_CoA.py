"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:84732 3-hydroxy fatty acyl-CoA

A 3-hydroxy fatty acyl-CoA is defined as a molecule that results from the formal condensation of
the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

Key features:
1. Contains a glycerol-3-phosphate moiety (G3P)
2. Contains a 3'-phosphoadenosine-5'-diphosphate moiety (ADP)
3. Contains a 3-hydroxy fatty acid chain attached via a thioester bond
4. The 3-hydroxy group is on the third carbon from the thioester
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for G3P moiety
    g3p_pattern = Chem.MolFromSmarts("OP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)O1")
    if not mol.HasSubstructMatch(g3p_pattern):
        return False, "Missing glycerol-3-phosphate moiety"
    
    # Look for ADP moiety
    adp_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O[C@H]1[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(adp_pattern):
        return False, "Missing 3'-phosphoadenosine-5'-diphosphate moiety"
    
    # Look for thioester bond
    thioester_pattern = Chem.MolFromSmarts("CS(=O)CC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester bond"
    
    # Look for 3-hydroxy fatty acid chain
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)CCCC")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Missing 3-hydroxy fatty acid chain"
    
    # Check position of hydroxy group
    hydroxy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1]
    for hydroxy in hydroxy_atoms:
        if hydroxy.GetTotalDegreeConnected() == 1:
            continue  # Terminal hydroxyl, not of interest
        atom3 = hydroxy.GetNeighbors()[0].GetNeighbors()[0]
        if atom3.GetAtomicNum() != 6:
            return False, "Hydroxy group not on 3rd carbon"
    
    # Count heavy atoms in fatty acid chain
    chain_atoms = [atom.GetIdx() for atom in mol.GetAtoms()
                   if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False]
    largest_chain = max(( len(Chem.SpanningTree(mol, chain_atoms, rootedAtomIdx=idx, isotopicAtoms=True)) for idx in chain_atoms ))

    if largest_chain < 6:
        return False, "Fatty acid chain too short"
    
    return True, "Contains 3-hydroxy fatty acid chain attached to CoA via thioester"