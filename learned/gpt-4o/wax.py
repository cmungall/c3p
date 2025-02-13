"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is considered a wax based on its SMILES string.
    A wax is characterized by long-chain molecules, typically containing esters,
    and being malleable at ambient temperatures.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a wax, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Check for long carbon chains (usually more than 14 carbons in a single chain)
    c_counts = [0] * (mol.GetNumAtoms() + 1)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = atom.GetNeighbors()
            c_neighbors = sum(1 for n in neighbors if n.GetAtomicNum() == 6)
            if c_neighbors == 2 or c_neighbors == 3:
                c_counts[atom.GetIdx()] = 1 + max(c_counts[n.GetIdx()] for n in neighbors if n.GetAtomicNum() == 6)
    
    max_chain = max(c_counts)
    if max_chain < 14:
        return False, f"Longest carbon chain is {max_chain}, which is too short for a wax"
    
    # Check for sufficient flexibility (many rotatable bonds)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Insufficient rotatable bonds ({n_rotatable}), waxes should be flexible"
    
    return True, f"Structure meets criteria for wax: Contains ester linkage, long carbon chain ({max_chain} carbons), and sufficient rotatable bonds ({n_rotatable})"