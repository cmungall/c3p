"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:36080 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is an sn-glycerol-3-phosphate with a single acyl group
    attached at either position 1 or position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol-3-phosphate backbone pattern (C-C-C with 2 oxygens and 1 phosphate group)
    backbone_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X3P]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No glycerol-3-phosphate backbone found"
    
    # Look for single ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be fatty acid"

    # Check molecular weight - monoacyl-sn-glycerol 3-phosphates typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for monoacyl-sn-glycerol 3-phosphate"

    # Count carbons, oxygens, and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 10:
        return False, "Too few carbons for monoacyl-sn-glycerol 3-phosphate"
    if o_count != 7:
        return False, "Must have exactly 7 oxygens (1 ester group, 1 phosphate group)"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"

    return True, "Contains glycerol-3-phosphate backbone with 1 fatty acid chain attached via ester bond"