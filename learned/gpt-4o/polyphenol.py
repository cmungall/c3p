"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol contains 2 or more benzene rings, two of which are substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a benzene ring
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    
    # Find all benzene rings in the molecule
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    # Collect benzene rings with hydroxy groups
    rings_with_hydroxy = set()
    for benzene_ring in benzene_matches:
        # Set to track if hydroxy found on this ring
        found_hydroxy = False
        for atom_idx in benzene_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    if len(neighbor.GetNeighbors()) == 1:  # Check if it's a -OH group
                        found_hydroxy = True
                        break  # One hydroxy per ring is sufficient
            if found_hydroxy:
                rings_with_hydroxy.add(tuple(sorted(benzene_ring)))
                break
    
    # Ensure at least 2 benzene rings have an OH substitution
    if len(rings_with_hydroxy) >= 2:
        return True, "Contains 2 or more benzene rings each with an OH group"
    
    return False, f"Only {len(rings_with_hydroxy)} benzene rings have a hydroxy group"