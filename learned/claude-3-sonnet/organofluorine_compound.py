"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:33839 organofluorine compound

An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for fluorine atoms
    has_fluorine = any(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms())
    
    # Check for carbon-fluorine bonds
    has_c_f_bond = any(atom.GetAtomicNum() == 9 and
                       any(nbr.GetAtomicNum() == 6 for nbr in atom.GetNeighbors())
                       for atom in mol.GetAtoms())
    
    if has_c_f_bond:
        return True, "Contains at least one carbon-fluorine bond"
    elif has_fluorine:
        return False, "Contains fluorine but no carbon-fluorine bonds"
    else:
        return False, "Does not contain any fluorine atoms"