"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a fatty acyl-CoA with a fatty acid chain of 6-12 C atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA structure: Pattern recognition for the adenine and phosphate components
    coa_pattern = Chem.MolFromSmarts("O=C2N(C=CN2)C3=NC=NC4=C3N=CN=C4N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A component"
    
    # Check for thioester linkage: -C(=O)-S-
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for medium-chain fatty acid length (6-12 carbons)
    fatty_acid_chain = Chem.MolFromSmarts("[C:1]([C])[C][C][C][C][C]")
    matches = mol.GetSubstructMatches(fatty_acid_chain)
    if not any(6 <= len(match) <= 12 for match in matches):
        return False, f"Fatty acid chain not in medium range (6-12 carbons), found {len(matches[0])} carbons" if matches else "No valid fatty acid chain found"

    return True, "Contains medium-chain fatty acyl-CoA structure with appropriate length"