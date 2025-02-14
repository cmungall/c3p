"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:38273 ether lipid

An ether lipid is defined as a lipid similar in structure to a glycerolipid but in which 
one or more of the carbon atoms on glycerol is bonded to an alkyl chain via an ether linkage,
as opposed to the usual ester linkage.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ether oxygen (-O-)
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Ether lipids must have at least one ether linkage
    if not ether_matches:
        return False, "No ether linkages found"
    
    # Ether lipids can have ester linkages, but not exclusively
    if ester_matches and not ether_matches:
        return False, "Found only ester linkages, need at least one ether linkage"
    
    # Check for long alkyl chains
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if not alkyl_chain_matches:
        return False, "No alkyl chains found"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short for ether lipid"
    
    return True, "Contains glycerol backbone with at least one ether linkage and long alkyl chains"