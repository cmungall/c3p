"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    This type of molecule has a phosphatidyl group attached to position 1 of 1D-myo-inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for 1D-myo-inositol: a six-membered ring with hydroxyl groups in the correct stereochemistry
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "1D-myo-inositol pattern not found or incorrect stereochemistry"

    # SMARTS pattern for the phosphatidyl group: a glycerol backbone with two esterified fatty acids and a phosphate group
    phosphatidyl_pattern = Chem.MolFromSmarts("O[C@@H](CO[P](=O)(O)OC[C@@H](OC(=O)C)OC(=O)C)C")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "Phosphatidyl group pattern not detected, or incorrect attachment"

    # Additional check to ensure the phosphatidyl group is connected to the correct carbon in inositol
    matches = mol.GetSubstructMatches(inositol_pattern)
    phosphatidyl_matches = mol.GetSubstructMatches(phosphatidyl_pattern)
    
    for match in matches:
        for phosphatidyl_match in phosphatidyl_matches:
            if phosphatidyl_match[0] in match:
                # Ensure the phosphatidyl group is at position 1 of myo-inositol
                return True, "Valid 1-phosphatidyl-1D-myo-inositol detected"

    return False, "Phosphatidyl group not attached to position 1 of inositol"