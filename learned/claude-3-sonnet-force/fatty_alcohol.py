"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:26372 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count alcohol groups (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) != 1:
        return False, f"Found {len(alcohol_matches)} alcohol groups, need exactly 1"
    
    # Check for aliphatic chain (only C, H, and O atoms)
    allowed_atoms = [6, 1, 8]  # C, H, O
    if any(atom.GetAtomicNum() not in allowed_atoms for atom in mol.GetAtoms()):
        return False, "Contains atoms other than C, H, and O"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 27:
        return False, f"Found {c_count} carbon atoms, need between 3 and 27"
    
    # Check for branching or unsaturation
    if mol.GetNumRotatableBonds() < c_count - 2:
        return True, "Unbranched and saturated fatty alcohol"
    elif mol.GetNumRotatableBonds() == c_count - 2:
        return True, "Unbranched and unsaturated fatty alcohol"
    else:
        return True, "Branched fatty alcohol"