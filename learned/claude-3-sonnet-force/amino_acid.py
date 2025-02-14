"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33709 amino acid
A carboxylic acid containing one or more amino groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is a carboxylic acid containing one or more amino groups.

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
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for primary amino group
    amino_pattern = Chem.MolFromSmarts("[NX3H2]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No primary amino group found"
    
    # Check for secondary and tertiary amino groups
    secondary_tertiary_amino_pattern = Chem.MolFromSmarts("[NX3H1,NX3H0]")
    secondary_tertiary_amino_matches = mol.GetSubstructMatches(secondary_tertiary_amino_pattern)
    
    # Count the total number of amino groups
    amino_group_count = len(mol.GetSubstructMatches(amino_pattern)) + len(secondary_tertiary_amino_matches)
    
    if amino_group_count < 1:
        return False, "No amino groups found"
    
    # Check for additional functional groups
    additional_groups = ["[SX2H0]", "[OX2H1]", "[OX1H0]", "[OX1H1]"]
    has_additional_groups = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in additional_groups)
    
    if has_additional_groups:
        return True, "Contains carboxylic acid and amino groups, with additional functional groups"
    else:
        return True, "Contains carboxylic acid and amino groups"