"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:36357 primary amine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound formally derived from ammonia by replacing
    one hydrogen atom by a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    
    # Check if any nitrogen atom has exactly one substituent (excluding H)
    for nitrogen in nitrogen_atoms:
        substituents = [neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetAtomicNum() != 1]
        if len(substituents) == 1:
            return True, "Molecule contains a nitrogen atom with one substituent (primary amine)"
    
    return False, "Molecule does not contain a primary amine group"