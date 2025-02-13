"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:38119 tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid (myristic acid)
with a hydroxy, amino, or quaternary ammonium group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for tetradecanoate ester moiety
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O[C,N]")

    # Check for tetradecanoate ester substructure
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester moiety found"

    # Check for glycerol backbone (glycerides)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    is_glyceride = mol.HasSubstructMatch(glycerol_pattern)

    # Check for phosphate group (phospholipids)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]")
    is_phospholipid = mol.HasSubstructMatch(phosphate_pattern)

    # Check for quaternary ammonium group (phospholipids)
    quat_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")
    has_quat_ammonium = mol.HasSubstructMatch(quat_ammonium_pattern)

    # Check for carnitine group
    carnitine_pattern = Chem.MolFromSmarts("C[N+](C)(C)CC([O-])=O")
    is_carnitine = mol.HasSubstructMatch(carnitine_pattern)

    # Classify based on structural patterns
    if is_glyceride:
        return True, "Contains glycerol backbone with tetradecanoate ester group"
    elif is_phospholipid or has_quat_ammonium:
        return True, "Contains phosphate or quaternary ammonium group with tetradecanoate ester group"
    elif is_carnitine:
        return True, "Contains carnitine group with tetradecanoate ester group"
    else:
        return False, "Does not match known structural patterns for tetradecanoate esters"