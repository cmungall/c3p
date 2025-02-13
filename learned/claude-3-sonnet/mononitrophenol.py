"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:51088 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find aromatic rings
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("c1ccccc1"))
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # Look for phenol/phenolate groups
    phenol_patterns = [
        Chem.MolFromSmarts("c1ccccc1[OH]"),    # phenol
        Chem.MolFromSmarts("c1ccccc1[O-]")     # phenolate
    ]
    
    valid_rings = []
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        oh_count = 0
        
        # Count OH/O- groups on this ring
        for pattern in phenol_patterns:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if any(idx in ring_atoms for idx in match):
                    oh_count += 1
        
        # Only consider rings with exactly one OH/O-
        if oh_count == 1:
            valid_rings.append(ring)
    
    if not valid_rings:
        return False, "No valid phenol group found (must have exactly one OH/O- group)"

    # Check for nitro groups
    nitro_patterns = [
        Chem.MolFromSmarts("[N+](=[O-])=O"),  # Charged form
        Chem.MolFromSmarts("[N+]([O-])=O"),    # Alternative charged form
        Chem.MolFromSmarts("N(=O)=O")          # Uncharged form
    ]
    
    nitro_matches = []
    for pattern in nitro_patterns:
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            nitro_matches.extend(matches)
    
    if not nitro_matches:
        return False, "No nitro group found"
    
    # Count unique nitro groups (some patterns might match the same group)
    unique_nitro_n = set(match[0] for match in nitro_matches)
    if len(unique_nitro_n) > 1:
        return False, f"Found {len(unique_nitro_n)} nitro groups, must have exactly one"

    # Check if the nitro group is attached to a valid phenol ring
    for ring in valid_rings:
        ring_atoms = set(ring)
        for nitro_n in unique_nitro_n:
            nitro_atom = mol.GetAtomWithIdx(nitro_n)
            neighbors = nitro_atom.GetNeighbors()
            
            for neighbor in neighbors:
                if neighbor.GetIdx() in ring_atoms:
                    return True, "Contains phenol ring with single nitro substituent"

    return False, "Nitro group not properly attached to phenol ring"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51088',
        'name': 'mononitrophenol',
        'definition': 'A nitrophenol that is phenol carrying a single nitro substituent at unspecified position.',
        'parents': ['CHEBI:33622']
    }
}