"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:33854 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is a primary alcohol where the alcoholic hydroxy group is attached
    to a carbon which is directly bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for aromatic primary alcohol
    # [#6;H2][OH] - primary alcohol carbon (C with two hydrogens) attached to OH
    # [$([#6;H2][OH])] - captures the primary alcohol group
    # [$(c)] - aromatic carbon
    # The pattern ensures the alcoholic carbon is connected to an aromatic carbon
    pattern = Chem.MolFromSmarts('[#6;H2][OH]')
    aromatic_primary_alcohol = False

    # Find all matches of the primary alcohol group
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No primary alcohol group found"

    for match in matches:
        alcohol_carbon_idx = match[0]  # Index of the alcoholic carbon
        alcohol_carbon = mol.GetAtomWithIdx(alcohol_carbon_idx)

        # Check if alcoholic carbon is connected to an aromatic carbon
        neighbors = alcohol_carbon.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIsAromatic():
                aromatic_primary_alcohol = True
                break  # Found aromatic attachment

        if aromatic_primary_alcohol:
            return True, "Molecule is an aromatic primary alcohol"

    if not aromatic_primary_alcohol:
        return False, "Primary alcohol is not connected to an aromatic ring"