"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are characterized by repeating isoprene units and a terminal alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for a single isoprene unit with variable connectivity
    isoprene_unit_pattern = Chem.MolFromSmarts("[C](=[C]([CH3])[*])[CH2][C]")
    matches = mol.GetSubstructMatches(isoprene_unit_pattern)
    num_isoprene_units = len(matches)

    if num_isoprene_units < 2:
      return False, f"Too few isoprene units, found {num_isoprene_units}"

    # Check for terminal alcohol (-C-OH) and ensure there is only one
    terminal_alcohol_pattern = Chem.MolFromSmarts("[!#1][OH]") # [!#1] is a non-H atom
    terminal_alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)

    if len(terminal_alcohol_matches) != 1 :
        return False, f"Incorrect number of terminal alcohol groups detected: {len(terminal_alcohol_matches)} (must be 1)"


    # Check for chain length (number of carbon atoms) - must be at least 5 * num_isoprene_units (roughly)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < (5 * num_isoprene_units):
        return False, "Too few carbons for the number of isoprene units."

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight is too low"


    return True, "Molecule has repeating isoprene units, a single terminal alcohol, sufficient chain length, and is compatible with polyprenol definition."