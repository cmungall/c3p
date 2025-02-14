"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(ncnc12)N2CC([C@H](O)C)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found"
    
    # Check for acyl group with thioester linkage
    acyl_thioester_pattern = Chem.MolFromSmarts("SCC(=O)C")
    if not mol.HasSubstructMatch(acyl_thioester_pattern):
        return False, "Acyl thioester linkage not found"
    
    # Check for 3-substitution at the propionyl group
    propionyl_substitution_pattern = Chem.MolFromSmarts("CC(C)C(=O)C[SH]")
    if not mol.HasSubstructMatch(propionyl_substitution_pattern):
        return False, "No 3-substituted propionyl group found"
    
    # Check for the presence of deprotonated phosphate groups (-4 charge)
    phosphate_deprotonation_pattern = Chem.MolFromSmarts("P([O-])(=O)(O)OP([O-])(=O)O")
    matches = mol.GetSubstructMatches(phosphate_deprotonation_pattern)
    if len(matches) < 2:
        return False, "Deprotonated phosphate groups indicating -4 charge not found"
    
    return True, "Molecule is a 3-substituted propionyl-CoA(4-)"