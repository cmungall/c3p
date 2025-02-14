"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: CHEBI:50974 isothiocyanate
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is an organosulfur compound with the general formula R-N=C=S.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isothiocyanate functional group patterns
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    has_isothiocyanate_group = mol.HasSubstructMatch(isothiocyanate_pattern)
    if not has_isothiocyanate_group:
        return False, "No isothiocyanate functional group found"
    
    # Check for carbon-nitrogen bond directly attached to isothiocyanate group
    direct_c_n_bond_pattern = Chem.MolFromSmarts("[C;X4]-[N;X2]=C=S")
    has_direct_c_n_bond = mol.HasSubstructMatch(direct_c_n_bond_pattern)
    
    if not has_direct_c_n_bond:
        return False, "No carbon-nitrogen bond directly attached to isothiocyanate group"
    
    return True, "Contains isothiocyanate functional group (R-N=C=S)"