"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: CHEBI:35831 cyclic fatty acid

A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("[C](=O)(O)")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No ring structure found"

    # Check for long carbon chain (fatty acid)
    longest_chain = AllChem.FindLongestChain(mol)
    if len(longest_chain.GetPaths()[0]) < 6:  # arbitrary threshold for fatty acid
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains a carboxylic acid group, ring structure, and long carbon chain"