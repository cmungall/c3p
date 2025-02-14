"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D is a seco-steroid with a characteristic broken B-ring, triene system,
    and hydroxyl groups at specific positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for characteristic steroid skeleton with broken B-ring
    # Note: this SMARTS pattern might be too broad, but it is a good starting point
    seco_steroid_pattern = Chem.MolFromSmarts("[C]1[CH]([CH2])[CH2][CH]([CH2])[CH2][CH2][CH](C)[CH2][CH2]1")
    if not mol.HasSubstructMatch(seco_steroid_pattern):
         return False, "Not a seco-steroid structure"

    # Check for triene system in the broken B-ring
    #  /C=C/C=C/C
    triene_pattern = Chem.MolFromSmarts("C=C-C=C-C")
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No triene system found"
    
    # Check for at least one hydroxyl group at C1, C3 or C25 (or combinations thereof)
    hydroxyl_pattern1 = Chem.MolFromSmarts("[CH0;R][OH1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern1)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl group found in characteristic positions"

    # Check for a side chain of at least 3 carbons
    side_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) < 1:
        return False, "No long enough side chain"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
         return False, f"Molecular weight out of range: {mol_wt}"

    return True, "Meets vitamin D structural criteria"