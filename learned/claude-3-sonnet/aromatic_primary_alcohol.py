"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:38092 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is defined as any primary alcohol in which the alcoholic
    hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.

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

    # Check for the presence of an aromatic ring
    has_aromatic_ring = mol.HasSubstructMatch(Chem.MolFromSmarts('c'))
    if not has_aromatic_ring:
        return False, "No aromatic ring found"

    # Check for the presence of a primary alcoholic group (-CH2-OH)
    has_primary_alcohol = mol.HasSubstructMatch(Chem.MolFromSmarts('[CH2,CH3]-[OH]'))
    if not has_primary_alcohol:
        return False, "No primary alcoholic group found"

    # Check if the alcoholic group is attached to a carbon bonded to an aromatic ring
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:  # Alcoholic oxygen
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 3:  # Primary carbon
                for c_neighbor in neighbor.GetNeighbors():
                    if c_neighbor.GetIsAromatic():
                        return True, "Contains an aromatic ring with a primary alcoholic group attached"

    return False, "Primary alcoholic group not attached to an aromatic ring"