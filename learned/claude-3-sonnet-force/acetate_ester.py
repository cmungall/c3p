"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35637 acetate ester
Any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetate group (-O-C(=O)C)
    acetate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CH3]")
    if not mol.HasSubstructMatch(acetate_pattern):
        return False, "No acetate group found"
    
    # Look for ester bond (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester bond found"
    
    # Check that the acetate group is part of an ester bond
    for match in ester_matches:
        acetate_match = mol.GetSubstructMatch(acetate_pattern)
        if any(atom_idx in match for atom_idx in acetate_match):
            return True, "Contains an acetate group as part of an ester bond"
    
    return False, "Acetate group not part of an ester bond"