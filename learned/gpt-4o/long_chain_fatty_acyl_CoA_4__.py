"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Long-chain fatty acyl-CoA(4-) are characteristic for having a long fatty acid chain
    attached to a Coenzyme A structure with a net charge of 4- due to deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A structure found"

    # Check for long-chain fatty acyl group (arbitrary definition, such as 12 or more carbons)
    # Long-chain feature represented as continuous chain of carbons, possibly interrupted by double bonds.
    fatty_acyl_pattern = Chem.MolFromSmarts("C[C,C,C=C,C=O]C(C)[C,C=C]")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    long_chain = any(len(match) >= 12 for match in fatty_acyl_matches)  # Simplified definition.
    if not long_chain:
        return False, "No long-chain fatty acyl group detected"

    # Ensure at least 4 phosphate groups with deprotonation (complex pattern of charge state)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-])([O-])O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 4:
        return False, "Insufficient phosphate groups (deprotonation assumed)"

    return True, "Identified as a long-chain fatty acyl-CoA(4-) structure"

# Example usage
smiles_example = "CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_long_chain_fatty_acyl_CoA_4__(smiles_example)
print(result, reason)