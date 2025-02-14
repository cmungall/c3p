"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:35569 secondary alcohol

A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has two other carbon atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary alcohol group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the number of hydroxyl groups
    num_hydroxyl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)

    # Check if there are any hydroxyl groups
    if num_hydroxyl == 0:
        return False, "No hydroxyl groups found"

    # SMARTS pattern for secondary alcohol group
    secondary_alcohol_pattern = Chem.MolFromSmarts("[CX4H2,CX4H1][OX2H]")

    # Check if the molecule has a substructure match for the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        # Check if the carbon atom is saturated and has exactly two other carbon neighbors
        for match in mol.GetSubstructMatches(secondary_alcohol_pattern):
            carbon_atom = mol.GetAtomWithIdx(match[0])
            if carbon_atom.GetHybridization() == Chem.HybridizationType.SP3 and len([n for n in carbon_atom.GetNeighbors() if n.GetAtomicNum() == 6]) == 2:
                return True, "Contains a secondary alcohol group"

    return False, "No secondary alcohol group found"