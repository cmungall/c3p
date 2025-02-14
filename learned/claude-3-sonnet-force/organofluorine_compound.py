"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:35480 organofluorine compound
An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Check for presence of fluorine atoms
    has_fluorine = any(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms())
    if not has_fluorine:
        return False, "No fluorine atoms present"
    
    # Check for carbon-fluorine bonds
    has_c_f_bond = any(bond.GetBeginAtom().GetAtomicNum() == 6 and
                       bond.GetEndAtom().GetAtomicNum() == 9
                       for bond in mol.GetBonds())
    
    if has_c_f_bond:
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bonds found"