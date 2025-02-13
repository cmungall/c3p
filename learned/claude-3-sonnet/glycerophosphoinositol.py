"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:18035 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with the polar alcohol inositol
    esterified to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone with phosphate group
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4][OX2][PX4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-phosphate backbone found"
    
    # Look for inositol group attached to phosphate
    inositol_pattern = Chem.MolFromSmarts("[PX4](-[OX2]1[CH2X4][CH1X4]([CH2X4][OH1X3])[OH1X3][CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])O1)-[OX2])")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol group attached to phosphate"
    
    # Look for 2 acyl/alkyl chains attached to glycerol backbone
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2X4]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 2:
        return False, f"Found {len(acyl_matches)} acyl/alkyl chains, need exactly 2"
    
    # Check molecular weight - glycerophosphoinositols typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycerophosphoinositol"
    
    return True, "Contains glycerol-phosphate backbone with inositol and 2 acyl/alkyl chains"