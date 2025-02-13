"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    A carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (C(=O)O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Look for amino group (N)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
        
    # Check the connectivity pattern to validate amino acid structure
    amino_acid_pattern = Chem.MolFromSmarts("[C;X4](N)(C(=O)O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No typical amino acid connectivity pattern found"

    return True, "Contains carboxylic acid group and amino group in the correct structure"