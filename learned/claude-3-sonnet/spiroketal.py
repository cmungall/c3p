"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:51633 spiroketal

A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ketal pattern (O-C-O)
    ketal_pattern = Chem.MolFromSmarts("[OX2]C([OX2])")
    if not mol.HasSubstructMatch(ketal_pattern):
        return False, "No ketal group found"
    
    # Find ring atoms
    ring_info = mol.GetRingInfo()
    
    # Check if ketal carbon is shared between two rings
    ketal_atoms = mol.GetSubstructMatches(ketal_pattern)
    ketal_carbon = ketal_atoms[0][1]
    ketal_rings = []
    for ring in ring_info.AtomRings():
        if ketal_carbon in ring:
            ketal_rings.append(ring)
    
    if len(ketal_rings) != 2:
        return False, "Ketal carbon not shared between two rings"
    
    # Check if ketal carbon is the only shared atom between the two rings
    ring1 = set(ketal_rings[0])
    ring2 = set(ketal_rings[1])
    shared_atoms = ring1.intersection(ring2)
    if len(shared_atoms) != 1 or list(shared_atoms)[0] != ketal_carbon:
        return False, "Ketal carbon not the only shared atom between rings"
    
    return True, "Contains a ketal group where the ketal carbon is shared between two rings"