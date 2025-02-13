"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a fundamental skeleton of four pyrrole nuclei connected by four methine groups forming a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for pyrrole
    pyrrole_pattern = Chem.MolFromSmarts("c1cc[nH]c1")
    
    # Check for exactly four pyrrole substructures
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Insufficient pyrrole units found: {len(pyrrole_matches)} (requires at least 4)"
    
    # SMARTS pattern for methine groups bridging
    methine_bridge_pattern = Chem.MolFromSmarts("C(c1cc[nH]c1)C(c1cc[nH]c1)C(c1cc[nH]c1)C")
    if not mol.HasSubstructMatch(methine_bridge_pattern):
        return False, "Methine bridge pattern not found"
    
    # Check for macrocyclic structure
    cycle_info = mol.GetRingInfo()
    if not cycle_info.IsAtomInRingOfSize(0, len(pyrrole_matches) + 4):  # Expect a large ring of size pyrroles + methines
        return False, "No large macrocyclic structure found"

    return True, "Contains four pyrrole units connected by methine bridges forming a macrocyclic structure"