"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    Fatty acyl-CoA is an acyl-CoA derived from the condensation of
    coenzyme A with a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the thioester linkage (C(=O)S) between acyl group and CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Enhanced SMARTS pattern for a fatty acyl chain, allowing for unsaturations
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)[C;!R][C;!R][C;!R][C;!R][C;!R]")  # At least a few carbons beyond the carbonyl
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No suitable fatty acyl chain found"

    # Simplified SMARTS pattern for coenzyme A structure, focusing on essential features
    coa_pattern = Chem.MolFromSmarts("C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n2cnc3c(N)ncnc23")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    return True, "Contains thioester linkage with fatty acyl chain and CoA structure"