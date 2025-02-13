"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[OX2]P([OX2H,OX1-])([OX2H,OX1-])=O")
    if len(mol.GetSubstructMatches(coa_pattern)) < 3:
        return False, "Missing phosphate groups characteristic of CoA"
        
    # Check for adenine
    adenine_pattern = Chem.MolFromSmarts("c1nc2c(n1)c(ncn2)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
        
    # Check for thioester bond
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester bond"
        
    # Count carbons in fatty acid chain
    # First find thioester carbon and traverse
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Cannot analyze fatty acid chain"
        
    # Count negative charges (should be 4)
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_count = len(mol.GetSubstructMatches(charge_pattern))
    if charge_count != 4:
        return False, f"Found {charge_count} negative charges, expected 4"
        
    # Analyze chain length
    # Use molecular descriptors to estimate chain length
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # CoA itself has about 23 carbons, so subtract to get fatty acid chain length
    chain_length = n_carbons - 23
    
    if chain_length < 4:
        return False, f"Fatty acid chain too short ({chain_length} carbons)"
    if chain_length > 14:
        return False, f"Fatty acid chain too long ({chain_length} carbons)"
        
    # Additional check for pantetheine arm
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm"

    return True, f"Medium-chain fatty acyl-CoA with approximately {chain_length} carbons in fatty acid chain"