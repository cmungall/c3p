"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for double bond between C5-C6
    # Note: The numbering in steroids starts from the A ring
    double_bond_5_6 = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]=,:[#6]~[#6]~3")
    if not mol.HasSubstructMatch(double_bond_5_6):
        return False, "No double bond between positions 5 and 6"

    # Check for 3-beta-hydroxyl group
    # We look for OH attached to the correct carbon in the A ring with correct stereochemistry
    beta_hydroxyl_3 = Chem.MolFromSmarts("[H][C@@H]1CC=C2C[C@@H](O)CC[C@]2")
    if not mol.HasSubstructMatch(beta_hydroxyl_3):
        return False, "No 3-beta-hydroxyl group found"

    # Additional validation - check for reasonable number of carbons (steroids typically have 17+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:
        return False, "Too few carbons for a steroid structure"

    return True, "Contains steroid core with 3-beta-hydroxyl group and double bond between C5-C6"