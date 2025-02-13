"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: carbonate ester
Definition: Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carbonate ester patterns
    # Pattern 1: R-O-C(=O)-O-R (symmetric or asymmetric)
    carbonate_pattern = Chem.MolFromSmarts("[O;X2][C;X3](=[O;X1])[O;X2]")
    
    # Pattern 2: R-O-C(=O)-OH (monocarbonate ester)
    monocarbonate_pattern = Chem.MolFromSmarts("[O;X2][C;X3](=[O;X1])[O;H1]")
    
    # Find matches
    carbonate_matches = mol.GetSubstructMatches(carbonate_pattern)
    monocarbonate_matches = mol.GetSubstructMatches(monocarbonate_pattern)
    
    if not (carbonate_matches or monocarbonate_matches):
        return False, "No carbonate ester group found"
    
    # For each carbonate match, verify it's not part of other functional groups
    for match in carbonate_matches:
        c_atom = match[1]  # The central carbon atom
        
        # Get neighboring atoms
        neighbors = [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(c_atom).GetNeighbors()]
        
        # Must have exactly 3 oxygen atoms connected
        if neighbors.count(8) != 3:  # 8 is atomic number for oxygen
            continue
            
        # Check if at least one oxygen is connected to carbon (organic group)
        organic = False
        for o_idx in [match[0], match[2]]:  # Check both O atoms
            o_atom = mol.GetAtomWithIdx(o_idx)
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom:
                    organic = True
                    break
            if organic:
                break
                
        if organic:
            return True, "Contains carbonate ester group (R-O-C(=O)-O-R or R-O-C(=O)-O-R')"
            
    # Check monocarbonate matches
    for match in monocarbonate_matches:
        c_atom = match[1]  # The central carbon atom
        
        # Get neighboring atoms
        neighbors = [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(c_atom).GetNeighbors()]
        
        # Must have exactly 3 oxygen atoms connected
        if neighbors.count(8) != 3:
            continue
            
        # Check if at least one oxygen is connected to carbon (organic group)
        for o_idx in [match[0]]:  # Check the non-OH oxygen
            o_atom = mol.GetAtomWithIdx(o_idx)
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom:
                    return True, "Contains monocarbonate ester group (R-O-C(=O)-OH)"
                    
    return False, "No valid carbonate ester group found"