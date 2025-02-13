"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: CHEBI:64426 beta-carbolines
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carboline(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines are pyridoindoles containing a beta-carboline skeleton and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyridoindole core
    pyridoindole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc1c3ccccc3n(c12)C")
    if not mol.HasSubstructMatch(pyridoindole_pattern):
        return False, "No pyridoindole core found"
    
    # Look for beta-carboline skeleton
    beta_carboline_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc1c3ccccc3n(c12)C4CNCC4")
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline skeleton"
    
    # Check for hydrogenated derivatives
    hydrogenated_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc1c3ccccc3n(c12)[C@H]4CCNCC4")
    if mol.HasSubstructMatch(hydrogenated_pattern):
        return True, "Hydrogenated beta-carboline derivative"
    
    return False, "Not a beta-carboline"