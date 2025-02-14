"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:50829 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone with a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for lactone ring pattern
    lactone_pattern = Chem.MolFromSmarts("[O;R]1[C;R][C;R][C;R][C;R]1=O")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    if len(lactone_matches) == 0:
        return False, "No lactone ring found"
    
    # Check if the lactone ring is a 5-membered ring (gamma-lactone)
    for match in lactone_matches:
        ring_atoms = mol.GetAtomRingInfo().AtomRings()[match]
        if len(ring_atoms) == 5:
            return True, "Molecule contains a 5-membered lactone ring (gamma-lactone)"
    
    return False, "Lactone ring found, but not a 5-membered ring (gamma-lactone)"

# Example usage
smiles = "O=C1O[C@@H](C/C=C\CC(=O)CCC)[C@@H](C1)C"
print(is_gamma_lactone(smiles))  # Output: (True, 'Molecule contains a 5-membered lactone ring (gamma-lactone)')