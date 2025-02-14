"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    This class can have a single acyl group at either position 1 or position 2 on the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for glycerol backbone attached with a phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(O)(O)=O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Pattern to match an ester linkage (-C(=O)O-) possibly on either sn-1 or sn-2 position
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[*:1]")
    ester_match_count = len(mol.GetSubstructMatches(ester_pattern))
    
    # Confirm the presence of only one ester linkage indicating one acyl group
    if ester_match_count != 1:
        return False, f"Incorrect number of ester linkages: Expected 1, found {ester_match_count}"
    
    return True, "Contains glycerol backbone with one acyl group and a phosphate group"