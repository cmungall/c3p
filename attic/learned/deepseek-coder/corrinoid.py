"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the corrin macrocycle pattern
    corrin_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]1")
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "No corrin macrocycle found"

    # Check for the presence of four pyrrole-like rings
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1cccc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole-like rings, need at least 4"

    # Check for the presence of three =C- groups and one direct carbon-carbon bond
    # This is a simplified check and may need refinement
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found {len(double_bond_matches)} =C- groups, need at least 3"

    # Check for the presence of a direct carbon-carbon bond linking alpha positions
    # This is a simplified check and may need refinement
    carbon_carbon_pattern = Chem.MolFromSmarts("[C]-[C]")
    carbon_carbon_matches = mol.GetSubstructMatches(carbon_carbon_pattern)
    if len(carbon_carbon_matches) < 1:
        return False, "No direct carbon-carbon bond found"

    return True, "Contains a corrin macrocycle with four pyrrole-like rings, three =C- groups, and one direct carbon-carbon bond"