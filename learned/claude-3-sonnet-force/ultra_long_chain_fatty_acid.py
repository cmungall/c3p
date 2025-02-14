"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:36977 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops
import numpy as np

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid (>C27) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count <= 27:
        return False, f"Chain length too short ({c_count} carbons)"
    
    # Check for fatty acid chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain found"
    
    # Check for substituents on the chain
    substituted_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[!C]")
    if mol.HasSubstructMatch(substituted_pattern):
        return False, "Carbon chain is substituted"
    
    # Count unsaturated carbons
    adjacency_matrix = rdmolops.GetAdjacencyMatrix(mol)
    unsaturated_bonds = np.sum(adjacency_matrix > 1)
    
    if unsaturated_bonds > 0:
        return True, f"Ultra-long-chain fatty acid with {c_count} carbons and {unsaturated_bonds} unsaturations"
    else:
        return True, f"Ultra-long-chain saturated fatty acid with {c_count} carbons"