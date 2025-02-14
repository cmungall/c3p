"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole core derived from L-tryptophan
    and a monoterpene unit derived from diisoprenoid building blocks connected together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole core
    indole_smarts = '[#6]-1:[#6]:[#6]:[#6]:[#6]:[#6]-1-[#7]-2:[#6]:[#6]:[#6]:[#6]:[#6]-2'  # simplified indole pattern
    indole_pattern = Chem.MolFromSmarts(indole_smarts)
    if indole_pattern is None:
        return False, "Invalid indole SMARTS pattern"

    # Check for indole core
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole core found"

    # Define SMARTS pattern for monoterpene unit (C10 terpenoid)
    monoterpene_smarts = '[#6]1[#6][#6][#6][#6][#6]1'  # simplified cyclohexene ring, may not be accurate
    monoterpene_pattern = Chem.MolFromSmarts(monoterpene_smarts)
    if monoterpene_pattern is None:
        return False, "Invalid monoterpene SMARTS pattern"

    # Attempt to find monoterpene units
    monoterpene_matches = mol.GetSubstructMatches(monoterpene_pattern)
    if len(monoterpene_matches) == 0:
        return False, "No monoterpene unit found"

    # Check if indole and monoterpene units are connected
    indole_atoms = mol.GetSubstructMatch(indole_pattern)
    monoterpene_atoms_list = mol.GetSubstructMatches(monoterpene_pattern)
    
    # Convert to sets for set operations
    indole_atom_set = set(indole_atoms)
    connected = False
    for monoterpene_atoms in monoterpene_atoms_list:
        monoterpene_atom_set = set(monoterpene_atoms)
        # Check for connectivity between indole and monoterpene units
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            if (begin_atom in indole_atom_set and end_atom in monoterpene_atom_set) or \
               (end_atom in indole_atom_set and begin_atom in monoterpene_atom_set):
                connected = True
                break
        if connected:
            break

    if not connected:
        return False, "Indole core and monoterpene unit are not connected"

    return True, "Contains an indole core connected to a monoterpene unit consistent with monoterpenoid indole alkaloids"