"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: CHEBI:36043 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is a benzene ring with exactly three chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"

    # Ensure there is only one benzene ring
    if len(benzene_matches) > 1:
        return False, "Multiple benzene rings found"

    # Check that all atoms are part of the benzene ring or its substituents
    benzene_atoms = set(benzene_matches[0])
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in benzene_atoms:
            # Check if the atom is a substituent of the benzene ring
            is_substituent = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in benzene_atoms:
                    is_substituent = True
                    break
            if not is_substituent:
                return False, "Molecule contains atoms not part of the benzene ring or its substituents"

    # Check for exactly three chlorine substituents on the benzene ring
    chlorine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Chlorine has atomic number 17
            # Check if the chlorine is attached to a benzene ring atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in benzene_atoms:
                    chlorine_count += 1
                    break

    if chlorine_count == 3:
        return True, "Contains a benzene ring with exactly three chlorine substituents"
    else:
        return False, f"Found {chlorine_count} chlorine substituents, need exactly 3"