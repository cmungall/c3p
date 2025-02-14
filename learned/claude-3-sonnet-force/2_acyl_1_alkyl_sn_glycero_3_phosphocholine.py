"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:18035 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule belongs to the class 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    based on its SMILES string. This class is defined as an alkyl,acyl-sn-glycero-3-phosphocholine
    with unspecified alkyl and acyl groups located at positions 1 and 2 respectively.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
    
    # Check for acyl and alkyl groups at positions 1 and 2
    acyl_pattern = Chem.MolFromSmarts("C(=O)[CH]")
    alkyl_pattern = Chem.MolFromSmarts("C[CH]")
    
    acyl_match = mol.GetSubstructMatches(acyl_pattern)
    alkyl_match = mol.GetSubstructMatches(alkyl_pattern)
    
    if len(acyl_match) != 1 or len(alkyl_match) != 1:
        return False, "Incorrect number of acyl/alkyl groups at positions 1 and 2"
    
    acyl_idx = acyl_match[0][1]
    alkyl_idx = alkyl_match[0][1]
    
    # Check if acyl and alkyl are at positions 1 and 2
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    
    if glycerol_atoms[0] != alkyl_idx and glycerol_atoms[2] != alkyl_idx:
        return False, "Alkyl group not at position 1"
    
    if glycerol_atoms[1] != acyl_idx:
        return False, "Acyl group not at position 2"
    
    return True, "Molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"