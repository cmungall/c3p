"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonic acid group (S(=O)(=O)-OH or S=O patterns)
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)[O-,OH]")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"
    
    # Check for carbon-sulfur bond to a long carbon chain or lipid
    c_s_bond_pattern = Chem.MolFromSmarts("[CX4,SX1]-[SX,OX2](=O)(=O)")
    if not mol.HasSubstructMatch(c_s_bond_pattern):
        return False, "No carbon-sulfur bond connecting sulfonic group to a lipid"

    return True, "Contains sulfonic acid residue joined by a carbon-sulfur bond to a lipid"