"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol structure with myo-configuration that has phosphate groups attached.

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

    # Check for cyclohexane core with exact hydroxyl position (stereochemistry included)
    # In myo-inositol, all hydroxyl groups are oriented in a specific manner
    myo_inositol_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1(O)")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "Molecule does not match the myo-inositol configuration"
    
    # Search for phosphate groups - note that phosphates can be complex with different charged states
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)[O-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"
    
    # Count number of phosphate groups, expecting at least one
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "Insufficient phosphate groups to qualify as myo-inositol phosphate"

    return True, "Molecule contains myo-inositol core with phosphate groups attached"