"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a macromolecule consisting of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts('[C;R]1-[C;R]-[C;R]-[C;R]-[C;R]-[O;R]-1')
    furanose_pattern = Chem.MolFromSmarts('[C;R]1-[C;R]-[C;R]-[C;R]-[O;R]-1')

    if pyranose_pattern is None or furanose_pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Find matches for pyranose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    # Find matches for furanose rings
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Total number of monosaccharide units
    total_units = len(pyranose_matches) + len(furanose_matches)

    if total_units > 10:
        return True, f"Contains {total_units} monosaccharide units linked glycosidically"
    else:
        return False, f"Contains only {total_units} monosaccharide units, need more than 10 for a polysaccharide"