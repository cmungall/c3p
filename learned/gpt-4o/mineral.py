"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved list of common metal and metalloid atomic numbers
    metal_and_metalloid_atomic_numbers = [
        3, 4, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 38, 39, 40, 50,
        51, 56, 58, 73, 74, 81, 82, 87, 88, 90 # Includes many commonly used metals and metalloids
    ]

    # Check for presence of metal ions or metalloids
    contains_metal_or_metalloid = any(atom.GetAtomicNum() in metal_and_metalloid_atomic_numbers for atom in mol.GetAtoms())
    if not contains_metal_or_metalloid:
        return False, "No metal ions or metalloids found in structure"
    
    # Check for presence of typical anions or groups in minerals
    mineral_group_smarts = [
        '[O-]S(=O)(=O)[O-]',  # sulfate
        'C(=O)([O-])[O-]',  # carbonate
        'P(=O)([O-])([O-])[O-]',  # phosphate
        '[N+]([O-])=O',  # nitrate
        'Cl',  # chloride
        '[O]',  # potential oxides/hydrates indicator
        '[S-]'  # sulfide
    ]
    
    for group_smarts in mineral_group_smarts:
        group_pattern = Chem.MolFromSmarts(group_smarts)
        if mol.HasSubstructMatch(group_pattern):
            # Distinguish from large organic species, e.g., checks such as few carbons or small structure true for ions
            carbon_atoms = [atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()]
            if sum(carbon_atoms) < 6:  # Small number of carbon atoms
                return True, "Contains metal ions/metalloids and suitable anion structure with low organic content."

    return False, "Does not match typical mineral structure"