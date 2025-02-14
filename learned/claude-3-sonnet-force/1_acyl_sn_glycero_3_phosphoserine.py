"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:17683 1-acyl-sn-glycero-3-phosphoserine

An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.

Examples:
- PS(16:1(9Z)/0:0): SMILES: P(OC[C@H](O)COC(=O)CCCCCCC/C=C\CCCCCC)(OC[C@H](N)C(O)=O)(O)=O
- PS(19:0/0:0): SMILES: P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
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
    
    # Look for phosphate group attached to the middle carbon of the glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][OX2]P([OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found on the glycerol backbone"
    
    # Look for serine head group
    serine_pattern = Chem.MolFromSmarts("[CH2][NX3H2]C([CX3](=[OX1]))[OX2H]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine head group found"
    
    # Look for acyl chain attached to the first carbon of the glycerol backbone
    acyl_pattern = Chem.MolFromSmarts("[CH2][CH2][OX2]C(=[OX1])[CX3]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl chain found"
    
    return True, "Contains glycerol backbone with phosphate, serine head group, and acyl chain"