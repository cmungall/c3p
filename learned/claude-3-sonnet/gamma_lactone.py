"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:35524 gamma-lactone
A lactone having a five-membered lactone ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone with a five-membered ring.

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
    
    # Look for lactone pattern (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("[O;R]1[C;R]=C[C;R][C;R]=O1")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    # Check if any of the matches have a 5-membered ring
    for match in lactone_matches:
        ring_atoms = mol.GetAtomRingInfo().AtomRings()[match]
        if len(ring_atoms) == 5:
            return True, "Contains a five-membered lactone ring"
    
    return False, "No five-membered lactone ring found"