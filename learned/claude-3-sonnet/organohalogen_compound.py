"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:33383 organohalogen compound
An organohalogen compound is defined as a compound containing at least one carbon-halogen bond 
(where X is a halogen atom).
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.

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
    
    # Check if molecule contains any halogen atoms
    halogen_atoms = ['Cl', 'Br', 'I', 'F']
    has_halogen = any(atom.GetSymbol() in halogen_atoms for atom in mol.GetAtoms())
    if not has_halogen:
        return False, "No halogen atoms present"
    
    # Check if any halogen atom is bonded to a carbon
    has_c_x_bond = any(atom.GetSymbol() in halogen_atoms and
                       any(nbr_atom.GetSymbol() == 'C' for nbr_atom in atom.GetNeighbors())
                       for atom in mol.GetAtoms())
    
    if has_c_x_bond:
        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bonds found"