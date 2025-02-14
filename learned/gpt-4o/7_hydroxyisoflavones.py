"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a specific isoflavone structure with a hydroxyl group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Isoflavone backbone pattern - aromatic rings with characteristic connectivity
    isoflavone_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)occ(-c1ccccc1)c2=O")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone backbone with 7-hydroxy group found"

    return True, "Molecule is a 7-hydroxyisoflavone with the correct structure"