"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is a steroid with a ketone group at position 3 and a double bond
    between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a 3-oxo group in steroids (carbonyl C=O at position 3)
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@H]")  # [C@H] is a chiral carbon indicator, could be adjusted if needed
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group (ketone) found on the third carbon"
    
    # SMARTS pattern for a Delta(1) double bond (C=C between first two carbons)
    double_bond_pattern = Chem.MolFromSmarts("C=C[C@H]")  # Includes chiral carbon subsequent to double bond
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No Delta(1) double bond (C=C between C1 and C2) found"
    
    # SMARTS pattern for a steroid core structure (four connected rings)
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC3C24C1")  # Simplified core structure
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have the steroid four-ring core structure"
    
    return True, "Matches 3-oxo-Delta(1) steroid structure"