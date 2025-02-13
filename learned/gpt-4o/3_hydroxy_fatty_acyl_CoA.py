"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is the result of condensation of the thiol group of coenzyme A
    with the carboxy group of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A pattern (simplified for matching the nucleotide portion and phosphates)
    coa_pattern = Chem.MolFromSmarts("O[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"

    # 3-hydroxy acyl pattern: C[CH](OH)C(=O) matching a 3-hydroxy group
    hydroxy_acyl_pattern = Chem.MolFromSmarts("C[CH](O)C(=O)S")
    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "3-hydroxy group with thioester linkage not found"

    return True, "Contains structure consistent with a 3-hydroxy fatty acyl-CoA"