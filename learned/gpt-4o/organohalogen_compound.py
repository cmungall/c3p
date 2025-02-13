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
    
    # Simplified SMARTS pattern for carbon-halogen bond
    organohalogen_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")  # Carbon bonded to any halogen
    
    # Check for matches using a general carbon-halogen bond pattern
    if mol.HasSubstructMatch(organohalogen_pattern):
        return True, "Contains carbon-halogen bond(s)"
    else:
        return False, "No carbon-halogen bond found"