"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: CHEBI:87851 aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alkylamino group
    alkylamino_pattern = Chem.MolFromSmarts("[N;X3]([C;H2,H3])[C;H2,H3]")
    alkylamino_matches = mol.GetSubstructMatches(alkylamino_pattern)
    if not alkylamino_matches:
        return False, "No alkylamino group found"

    # Check for aromatic ring directly attached to alkyl group
    aromatic_ring_pattern = Chem.MolFromSmarts("*@[ar]")
    for match in alkylamino_matches:
        n_atom = mol.GetAtomWithIdx(match[0])
        neighbors = n_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.IsInRing() and mol.GetAtomWithIdx(neighbor.GetIdx()).GetIsAromatic():
                return True, "Contains alkylamino group with aromatic ring directly attached"

    return False, "No aromatic ring directly attached to alkylamino group"