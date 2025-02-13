"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:27025 mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or structural functions in the human body.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound containing essential mineral elements.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is inorganic
    if any(atom.GetAtomicNum() in (6, 8, 1) for atom in mol.GetAtoms()):
        return False, "Contains organic atoms (C, H, O)"

    # Check for essential mineral elements
    mineral_elements = [3, 9, 11, 12, 15, 16, 17, 19, 20, 25, 26, 27, 28, 29, 30, 38, 48, 53, 56, 57]
    elements_present = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    mineral_elements_present = elements_present.intersection(mineral_elements)
    if not mineral_elements_present:
        return False, "Does not contain any essential mineral elements"

    # Check for common mineral nutrient anions
    anion_patterns = [
        Chem.MolFromSmarts("[O-]P([O-])([O-])=O"),  # phosphate
        Chem.MolFromSmarts("[O-]S([O-])(=O)=O"),    # sulfate
        Chem.MolFromSmarts("[Cl-]"),                # chloride
        Chem.MolFromSmarts("[F-]"),                 # fluoride
        Chem.MolFromSmarts("[O-]C([O-])=O"),        # carbonate
        Chem.MolFromSmarts("[O-][N+]([O-])=O"),     # nitrate
    ]
    anion_matches = [mol.HasSubstructMatch(pattern) for pattern in anion_patterns]
    if not any(anion_matches):
        return False, "Does not contain common mineral nutrient anions"

    return True, "Inorganic compound containing essential mineral elements and mineral nutrient anions"