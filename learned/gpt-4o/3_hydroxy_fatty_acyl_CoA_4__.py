"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

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
    
    # Look for the 3-hydroxy fatty acid structure
    # Acyl-CoA feature: Typically a long hydrocarbon chain with a hydroxyl group; examine an example SMILES pattern
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)C(=O)SC")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acyl structure found"
    
    # Look for Coenzyme A structure
    # Here, assume CoA includes the pantetheine unit and the nucleotide adenine component 
    coa_motif = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1O[C@H](n2cnc3c2ncnc3N)C1O")
    if not mol.HasSubstructMatch(coa_motif):
        return False, "CoA motif not present"
    
    # Verify deprotonation; assume negatively charged oxygen by phosphorus
    deprotonated_phosphate = Chem.MolFromSmarts("P([O-])(=O)([O-])")
    if len(mol.GetSubstructMatches(deprotonated_phosphate)) < 2:
        return False, "Not enough deprotonated phosphate groups found"
    
    return True, "Molecule classified as 3-hydroxy fatty acyl-CoA(4-)"