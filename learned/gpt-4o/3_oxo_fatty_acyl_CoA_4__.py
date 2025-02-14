"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class features a 3-oxo group, a fatty acyl chain, and a Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the Coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("OP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Check for the 3-oxo group pattern
    oxo_group_pattern = Chem.MolFromSmarts("CC(=O)")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"

    # Check for long fatty acyl chain
    acyl_chain_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]-[#6]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No long fatty acyl chain found"

    return True, "Valid 3-oxo-fatty acyl-CoA(4-) structure identified"

# Example usage:
smiles = "CCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles)
print(f"Result: {result}, Reason: {reason}")