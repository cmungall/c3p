"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determine if a molecule is an organometalloidal compound based on its SMILES string.
    Organometalloidal compounds have bonds between metalloids (e.g., As) and organic carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Patterns
    arsenic_carbon_pattern = Chem.MolFromSmarts("[As][C]")
    organoarsenic_pattern = Chem.MolFromSmarts("[As](=O)(O)([C])")  # Expand to cover varied functional arsenic configurations
    aromatic_arsenic_pattern = Chem.MolFromSmarts("c[As]")
    organometal_multiple_arsenic_pattern = Chem.MolFromSmarts("[C][As][As][C]")  # Multiple arsenic atoms involving carbon
    
    found_patterns = []

    # Check for each pattern and record reason
    if mol.HasSubstructMatch(arsenic_carbon_pattern):
        found_patterns.append("Contains arsenic-carbon bonds")
    if mol.HasSubstructMatch(organoarsenic_pattern):
        found_patterns.append("Contains functional organyl-arsenic component")
    if mol.HasSubstructMatch(aromatic_arsenic_pattern):
        found_patterns.append("Contains arsenic in an aromatic context")
    if mol.HasSubstructMatch(organometal_multiple_arsenic_pattern):
        found_patterns.append("Contains multiple arsenic atoms with organic linkages")

    if found_patterns:
        return True, " and ".join(found_patterns) + " indicative of organometalloidal compounds"

    return False, "Does not meet organometalloidal structural criteria"