"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group at the beta- or 3-position relative to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Position of carboxylic acid's carbon
    carb_ketone_pos = carboxylic_acid_matches[0][0]

    # Identify hydroxyl group attached to third carbon from carboxyl carbon
    for atom in mol.GetAtoms():
        if atom.GetIdx() == carb_ketone_pos + 2:  # Check the third carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Check if oxygen (OH group)
                    if any(nbr.GetAtomicNum() == 1 for nbr in neighbor.GetNeighbors()):  # Ensure it's an OH group
                        return True, "Contains hydroxyl group at the 3-position and carboxylic acid group"
    
    return False, "Hydroxyl group not found at the 3-position"