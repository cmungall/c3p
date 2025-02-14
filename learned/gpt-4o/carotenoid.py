"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 60:
        return False, f"Carotenoids typically have around 40 carbons, found {c_count}"

    # Look for a significant number of conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=C")
    conjugated_bonds = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_bonds) < 10:
        return False, "Insufficient conjugated double bonds"

    # Check for multiple oxygen atoms indicative of functional diversity (without strict need for OH)
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen < 1:
        return False, "No oxygen atoms found"

    # Assess for ring presence, allowing for linear structures or complexity
    ring_info = mol.GetRingInfo().BondRings()
    if len(ring_info) > 0 or c_count >= 40:
        return True, "Matches the structural criteria of carotenoids"

    return False, "No ring structures or insufficient complexity to classify as carotenoid"