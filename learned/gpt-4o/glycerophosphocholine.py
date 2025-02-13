"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine is characterized by a glycerol backbone, two fatty acids or ether-linked aliphatic side chains,
    and a phosphocholine moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (launches proper positions)
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("O[P@](=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check if there are two long chains (fatty acid esters or ethers)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("COC")
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))

    if ester_matches + ether_matches < 2:
        return False, f"Insufficient chain linkages found (esters: {ester_matches}, ethers: {ether_matches})"

    return True, "Contains glycerol backbone, phosphocholine group, and at least two aliphatic side chains as ester or ether linkages"