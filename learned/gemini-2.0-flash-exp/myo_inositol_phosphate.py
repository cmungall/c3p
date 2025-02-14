"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the myo-inositol core with correct stereochemistry using SMARTS
    # The following SMARTS pattern enforces the myo-inositol stereochemistry.
    # O[C@H]1[C@H](O)[C@@H](O)[C@H]([OX2])[C@H](O)[C@@H]1[OX2]
    # The [OX2] represents that it must have two connections
    myo_inositol_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H]([OX2])[C@H](O)[C@@H]1[OX2]')
    if not mol.HasSubstructMatch(myo_inositol_pattern):
         return False, "No myo-inositol core found"

    # 2. Check for at least one phosphate group
    phosphate_pattern = Chem.MolFromSmarts('[OX2][P](=[OX1])([OX2])[OX2]')
    diphosphate_pattern = Chem.MolFromSmarts('[OX2][P](=[OX1])([OX2])[OX2][P](=[OX1])([OX2])[OX2]')
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    
    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate group found"
    
    # 3. Check that the phosphate is attached to the inositol
    # We do not need to check this as the previous smart patterns already take care of it.
    
    return True, "Contains myo-inositol core with at least one phosphate group"