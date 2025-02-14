"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:51601 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound based on its SMILES string.
    A volatile organic compound is any organic compound having an initial boiling point
    less than or equal to 250 °C (482 °F) measured at a standard atmospheric pressure
    of 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a volatile organic compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for organic compound (contains carbon atoms)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        return False, "No carbon atoms found, not an organic compound"
    
    # Calculate molecular descriptors related to boiling point
    mol_weight = Descriptors.MolWt(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    topological_polar_surface_area = Descriptors.TPSA(mol)
    
    # Define a heuristic rule for volatility based on descriptors
    # This is a simplified example, and the rule can be improved or replaced with a machine learning model
    if mol_weight < 150 and num_rotatable_bonds < 10 and topological_polar_surface_area < 50:
        return True, f"Molecular descriptors suggest a volatile organic compound: MW={mol_weight:.2f}, RotB={num_rotatable_bonds}, TPSA={topological_polar_surface_area:.2f}"
    else:
        return False, f"Molecular descriptors suggest a non-volatile organic compound: MW={mol_weight:.2f}, RotB={num_rotatable_bonds}, TPSA={topological_polar_surface_area:.2f}"