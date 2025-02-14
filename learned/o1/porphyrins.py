"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is characterized by a macrocyclic structure composed of four pyrrole rings connected via methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if the ring is 16-membered
        if len(ring) == 16:
            # Count the number of nitrogen atoms in the ring
            num_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            # Check if there are exactly four nitrogen atoms
            if num_nitrogen == 4:
                return True, "Contains porphyrin core: 16-membered ring with 4 nitrogen atoms"
    
    return False, "Does not contain porphyrin core macrocyclic structure"