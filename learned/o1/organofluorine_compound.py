"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is a compound containing at least one carbon-fluorine bond.

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
    
    # Initialize flag
    has_c_f_bond = False
    
    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        atomic_nums = {atom1.GetAtomicNum(), atom2.GetAtomicNum()}
        
        # Check if bond is between carbon (6) and fluorine (9)
        if atomic_nums == {6, 9}:
            has_c_f_bond = True
            break  # No need to check further
    
    if has_c_f_bond:
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "Does not contain any carbon-fluorine bonds"