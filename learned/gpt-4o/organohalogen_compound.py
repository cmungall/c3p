"""
Classifies: CHEBI:17792 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more specific SMARTS patterns for carbon-halogen bonds
    fragrance_halogen_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")  # Specifically ensuring carbon is singly bonded to halogen
    alkyl_halide_pattern = Chem.MolFromSmarts("[CX4]-[F,Cl,Br,I]")      # Alkyl halide - saturated carbon
    aryl_halide_pattern = Chem.MolFromSmarts("[cX3]-[F,Cl,Br,I]")      # Aryl halide - aromatic carbon
    
    # Check for matches
    if (mol.HasSubstructMatch(fragrance_halogen_pattern) or
        mol.HasSubstructMatch(alkyl_halide_pattern) or
        mol.HasSubstructMatch(aryl_halide_pattern)):
        return True, "Contains carbon-halogen bond(s)"
    else:
        return False, "No carbon-halogen bond found"