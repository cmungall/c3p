"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group found"

    # Look for aromatic pattern
    aromatic = False
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic = True
            break

    # Check that nitrile group is not in an aromatic context
    # and is part of an aliphatic chain by ensuring absence of aromatic bonds in nitrile carbon's neighbors
    for match in mol.GetSubstructMatches(nitrile_pattern):
        carbon_idx, nitrogen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # If the nitrile carbon itself is marked aromatic, it's not aliphatic
        if carbon_atom.GetIsAromatic():
            return False, "Nitrile group is part of an aromatic structure"

        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return False, "Nitrile group is attached to an aromatic ring"

    return True, "Nitrile group is derived from an aliphatic compound"