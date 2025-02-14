"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    A short-chain fatty acyl-CoA is a fatty acyl-CoA formed by the condensation
    of the thiol group of coenzyme A with a short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the Coenzyme A core structure pattern (simplified)
    coenzymeA_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCS")
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "No coenzyme A moiety found"
    
    # Define the thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for short-chain fatty acid component (2-6 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("C(CCCC)C")
    matches = mol.GetSubstructMatches(carbon_chain_pattern)
    
    # Count the total number of carbons in the chain
    if not any(2 <= sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6) <= 6 for match in matches):
        return False, "No short-chain fatty acid component found"
    
    # If all checks pass, classify as short-chain fatty acyl-CoA
    return True, "Contains a CoA moiety joined by thioester bond to a short-chain fatty acid"