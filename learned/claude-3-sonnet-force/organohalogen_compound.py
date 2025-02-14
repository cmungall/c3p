"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:35922 organohalogen compound
A compound containing at least one carbon-halogen bond (where X is a halogen atom).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound is defined as a compound containing at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the list of atoms
    atoms = mol.GetAtoms()
    
    # Check if any atom is carbon bonded to a halogen
    for atom in atoms:
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
                    return True, "Contains at least one carbon-halogen bond"
    
    # No carbon-halogen bond found
    return False, "No carbon-halogen bond found"