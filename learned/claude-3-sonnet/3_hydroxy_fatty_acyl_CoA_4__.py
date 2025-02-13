"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:15529 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "No adenine moiety found"

    # Check for 4 negative charges (phosphate groups)
    phosphate = Chem.MolFromSmarts("[O-]P([O-])(=O)O")
    if len(mol.GetSubstructMatches(phosphate)) < 2:
        return False, "Missing required phosphate groups"

    # Check for thioester linkage
    thioester = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy group on fatty acid chain
    hydroxy_pattern = Chem.MolFromSmarts("CC(O)CC(=O)S")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group found in correct position"

    # Check for pantetheine part
    pantetheine = Chem.MolFromSmarts("CC(C)(CO)[CH]OC(=O)CCNC(=O)CCNC")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Check for ribose
    ribose = Chem.MolFromSmarts("C1C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(ribose):
        return False, "Missing ribose moiety"

    # Count carbons in fatty acid chain (should be at least 4)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)SCCNC")
    matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not matches:
        return False, "No fatty acid chain found"

    return True, "Contains CoA(4-) moiety with 3-hydroxy fatty acid chain attached via thioester bond"