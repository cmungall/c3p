"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: CHEBI:51205 beta-carboline
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_carboline(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.

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
    
    # Look for pyridoindole core (fused pyridine and indole rings)
    pyridoindole_pattern = Chem.MolFromSmarts("[*]1nc2ccccc2c3[nH]ccc13")
    if not mol.HasSubstructMatch(pyridoindole_pattern):
        return False, "No pyridoindole core found"
    
    # Look for beta-carboline skeleton (fused pyridine, indole, and cyclohexene rings)
    beta_carboline_pattern = Chem.MolFromSmarts("[*]1nc2ccccc2c3c1[nH]c4ccccc34")
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline skeleton"
    
    # Check for hydrogenated derivatives
    hydrogenated_pattern = Chem.MolFromSmarts("[*]1nc2ccccc2c3c1[nH]c4ccccc34")
    if mol.HasSubstructMatch(hydrogenated_pattern):
        return True, "Hydrogenated derivative of beta-carboline"
    
    return False, "Not a beta-carboline"