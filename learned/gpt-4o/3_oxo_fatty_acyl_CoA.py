"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is a fatty acyl-CoA with an oxo group on the third carbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the CoA pattern
    coa_pattern = Chem.MolFromSmarts("[P](=O)([O-])[O-].[C@H](O)C(C)(C)COP(=O)(O)O[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Look for 3-oxo linkage pattern
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts("CC(=O)C(=O)S[C,N]") # Pattern to identify a 3-oxo-fatty acyl
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "No 3-oxo-fatty acyl linkage found"

    return True, "Contains CoA structure and 3-oxo-fatty acyl linkage"

# Example usage
example_smiles = "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
print(is_3_oxo_fatty_acyl_CoA(example_smiles))