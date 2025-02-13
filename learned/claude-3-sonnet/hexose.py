"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18107 hexose
A hexose is any six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count < 5 or o_count > 6:
        return False, "Incorrect number of carbon (6) or oxygen (5-6) atoms"

    # Check for aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("[CH2](C=O)")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    if len(aldehyde_matches) == 1:
        return True, "Molecule contains an aldehyde group at position 1 (aldohexose)"
    elif len(ketone_matches) == 1:
        return True, "Molecule contains a ketone group at position 2 (ketohexose)"
    else:
        return False, "No aldehyde or ketone group found in the correct position"

    return False, "Molecule does not meet hexose criteria"