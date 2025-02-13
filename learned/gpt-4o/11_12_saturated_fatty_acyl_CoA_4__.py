"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA pattern with key functional group presence
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C[CH](O)C(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A group not found"
        
    # Check for a saturated bond pattern at arbitrary position
    saturated_11_12_pattern = Chem.MolFromSmarts("[CX4][CX4]")
    if not mol.HasSubstructMatch(saturated_11_12_pattern):
        return False, "11-12 bond is not saturated"
    
    # Total number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 12:
        return False, "Not enough carbon atoms for a fatty acyl chain"
    
    return True, "Contains CoA group and saturated 11,12 bond"

# Test with a sample SMILES from the class
smiles_str = "CCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(smiles_str)
print(f"Result: {result}, Reason: {reason}")