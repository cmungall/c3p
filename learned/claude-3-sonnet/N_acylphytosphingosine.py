"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:84670 N-acylphytosphingosine

A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phytosphingosine backbone pattern
    phytosphingosine_pattern = Chem.MolFromSmarts("[NH1][CX4]([CH2][CH2][CH2][CH2][CH2])[CX4]([CH1]([OH1])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2])[CH2][OH1]")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "Phytosphingosine backbone not found"
    
    # Look for acyl chain attached to N atom
    acyl_pattern = Chem.MolFromSmarts("[NH1][CX3](=[OX1])[CX3]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to nitrogen"
    
    # Check for long acyl chain (>= 6 carbons)
    long_acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]")
    if not mol.HasSubstructMatch(long_acyl_pattern):
        return False, "Acyl chain too short (< 6 carbons)"
    
    return True, "Contains phytosphingosine backbone with acyl group attached to nitrogen"