"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol molecule where two of the three hydroxyl groups
    are substituted with acyl (ester-linked), alkyl (ether-linked), or alk-1-enyl (vinyl ether-linked) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        substituted_positions = 0  # Reset for each glycerol backbone

        # For each carbon in the glycerol backbone
        for c_idx in [c1_idx, c2_idx, c3_idx]:
            carbon = mol.GetAtomWithIdx(c_idx)
            # Find oxygens attached to this carbon
            oxygen_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if not oxygen_neighbors:
                continue  # No oxygen attached to this carbon (unlikely in glycerol)
            oxygen = oxygen_neighbors[0]
            # Check if oxygen is substituted
            heavy_atom_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetAtomicNum() != 1]
            if len(heavy_atom_neighbors) > 1:
                # Oxygen connected to another heavy atom besides glycerol carbon (substituted)
                substituted_positions += 1

        if substituted_positions == 2:
            return True, "Molecule is a diradylglycerol with two substituted positions"

    return False, f"Found {substituted_positions} substituted positions, expected 2"