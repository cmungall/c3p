"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:18311 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile has nitrile groups (-C#N) attached to aliphatic carbons (non-aromatic).

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

    # Find all nitrile groups (C#N)
    nitrile_pattern = Chem.MolFromSmarts('[C]#[N]')
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not nitrile_matches:
        return False, "No nitrile groups found"
    
    # Check each nitrile's carbon atom for aromaticity
    for match in nitrile_matches:
        # Match tuple contains (carbon_idx, nitrogen_idx)
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetIsAromatic():
            return False, f"Nitrile group attached to aromatic carbon at position {match[0]+1}"

    return True, "All nitrile groups are attached to aliphatic carbons"