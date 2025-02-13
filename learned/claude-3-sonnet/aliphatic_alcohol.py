"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:15892 aliphatic alcohol

An aliphatic alcohol is an alcohol derived from an aliphatic compound.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Check for aliphatic backbone (no cycles, no aromatic rings)
    aliphatic_pattern = Chem.MolFromSmarts("[!#1!#6!r]")  # Not H, not C, not in ring
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if aliphatic_matches:
        return False, "Contains non-aliphatic atoms or rings"

    # Check for presence of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 1:
        return False, "No carbon atoms found"

    # Count rotatable bonds to check for long aliphatic chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Insufficient rotatable bonds for aliphatic chains"

    return True, "Contains an alcohol group and aliphatic backbone"