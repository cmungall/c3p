"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.

    A mineral nutrient is a mineral that is an inorganic nutrient which must be ingested 
    and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or 
    structural functions in the human body.

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

    # List of essential mineral elements (expanded)
    essential_minerals = {
        "Ca",  # Calcium
        "P",   # Phosphorus
        "K",   # Potassium
        "Na",  # Sodium
        "Cl",  # Chlorine
        "Mg",  # Magnesium
        "S",   # Sulfur
        "Fe",  # Iron
        "Mn",  # Manganese
        "Cu",  # Copper
        "I",   # Iodine
        "Zn",  # Zinc
        "Se",  # Selenium
        "Mo",  # Molybdenum
        "Cr",  # Chromium
        "Co",  # Cobalt
        "F",   # Fluorine
        "Ba",  # Barium
        "Sb",  # Antimony
        "La",  # Lanthanum
        "Cs",  # Cesium
        "Pd",  # Palladium
        "Al",  # Aluminium
        "Si"   # Silicon
    }

    # Flag to indicate presence of essential mineral not bonded to carbon
    has_essential_mineral = False

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in essential_minerals:
            is_bonded_to_carbon = False
            # Check if the essential mineral is directly bonded to carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "C":
                    is_bonded_to_carbon = True
                    break
            if not is_bonded_to_carbon:
                has_essential_mineral = True
                # No need to check further if we found an essential mineral
                break

    if not has_essential_mineral:
        return False, "Does not contain essential mineral elements or elements are bonded to carbon"

    # Additional checks to exclude complex organic molecules
    # Check if the molecule contains rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains ring structures; likely not a simple mineral nutrient"

    # Count total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 20:
        return False, f"Contains {c_count} carbon atoms; too many for a simple mineral nutrient"

    return True, "Contains essential mineral elements not bonded to carbon; no rings; limited carbon atoms"