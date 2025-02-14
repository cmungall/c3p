"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True if molecule is an amino acid, False otherwise and the reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for amino and carboxylic acid groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]") # includes primary, secondary and tertiary amines
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1]") # covers carboxylate and carboxylic acid forms

    # Check for presence of amino group
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Check for presence of carboxylic acid group
    if not mol.HasSubstructMatch(carboxyl_pattern):
          return False, "No carboxylic acid group found"

    # Check if it's an alpha amino acid by looking for 1 carbon distance between the carboxyl carbon and the amino nitrogen
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4][CX3](=[OX1])[OX2]")

    #if an alpha amino acid is found, then return True and the appropriate reason.
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Contains both amino and carboxylic acid groups with alpha carbon"
    
    # if not an alpha amino acid return true because it contains an amino and carboxylic acid
    return True, "Contains both amino and carboxylic acid groups"