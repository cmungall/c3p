"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is an oxo steroid that has an oxo (=O) substituent at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid nucleus (four fused rings)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C(C1)CC4CC(C(C2)C3)C4')  # Basic steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Molecule does not have a steroid nucleus"

    # Identify the rings in the molecule
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:
        return False, "Molecule does not have the required four rings of a steroid nucleus"

    # Separate six-membered and five-membered rings
    six_membered_rings = [ring for ring in ssr if len(ring) == 6]
    five_membered_rings = [ring for ring in ssr if len(ring) == 5]

    if len(six_membered_rings) < 3 or len(five_membered_rings) < 1:
        return False, "Steroid nucleus must have three six-membered rings and one five-membered ring"

    # Define SMARTS pattern for ketone (=O) group
    ketone_pattern = Chem.MolFromSmarts('C=O')
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone (=O) group found in molecule"

    # Check if any ketone is in ring C (the third six-membered ring)
    ring_C = six_membered_rings[2]  # Assuming ring C is the third six-membered ring
    ketone_in_ringC = False
    for match in ketone_matches:
        carbon_idx = match[0]
        if carbon_idx in ring_C:
            ketone_in_ringC = True
            break

    if not ketone_in_ringC:
        return False, "Ketone group is not located on ring C (position 11)"

    return True, "Molecule has a steroid nucleus with a ketone group at position 11"