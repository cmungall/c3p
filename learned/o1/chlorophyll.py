"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll

Definition: A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is characterized by a magnesium chlorin core with a fifth ring and side chains,
    often including a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens (may help with valence issues)
    mol = Chem.AddHs(mol)

    # Check for magnesium atom (atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"

    # For each magnesium atom
    for mg in mg_atoms:
        # Find nitrogen atoms bonded to magnesium
        neighbors = mg.GetNeighbors()
        nitrogen_neighbors = [atom for atom in neighbors if atom.GetAtomicNum() == 7]
        if len(nitrogen_neighbors) >= 4:
            # Check if nitrogen atoms are part of rings
            ring_nitrogens = [n for n in nitrogen_neighbors if n.IsInRing()]
            if len(ring_nitrogens) >= 4:
                # Potential chlorophyll
                return True, "Contains magnesium coordinated to 4 ring nitrogen atoms"
    return False, "Magnesium not coordinated to 4 ring nitrogen atoms"

__metadata__ = {
    'chemical_class': {
        'name': 'chlorophyll',
        'definition': 'A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.'
    }
}