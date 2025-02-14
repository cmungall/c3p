"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:XXXXXX 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a steroid with a hydroxy group at position 16 in the beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid backbone SMARTS pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C)(C4CCC3C2C1)')  # Simplified pattern
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define SMARTS pattern for 16beta-hydroxy group
    # Position 16 beta-hydroxy group on steroid backbone
    # This pattern may need adjustment based on exact stereochemistry
    hydroxy16beta_pattern = Chem.MolFromSmarts('[C@H](O)[C@@]1(CC[C@H]2[C@H](C)CC[C@@]3(C)C(=CC=C4)C4CCC3C2C1)C')  # Simplified

    if mol.HasSubstructMatch(hydroxy16beta_pattern):
        return True, "Steroid backbone with 16beta-hydroxy group found"

    # Attempt to find the atom corresponding to position 16
    # For accuracy, we might need to manually map the atoms based on known examples
    # Here, we'll search for a chiral center with hydroxyl group in beta position

    # Find all chiral carbons with hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.HasProp('_ChiralityPossible'):
            neighbors = atom.GetNeighbors()
            has_oxygen = any(neighbor.GetAtomicNum() == 8 for neighbor in neighbors)
            if has_oxygen:
                chiral_tag = atom.GetChiralTag()
                if chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    # Assuming CHI_TETRAHEDRAL_CCW corresponds to beta-configuration
                    return True, "Found chiral carbon with beta-hydroxy group, possible 16beta-hydroxy steroid"

    return False, "No 16beta-hydroxy group with beta-configuration found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXXX',
        'name': '16beta-hydroxy steroid',
        'definition': 'A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here
    }
}