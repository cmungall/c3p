"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: fatty acid anion
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one carboxylate group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Molecule has {len(carboxylate_matches)} carboxylate groups, expected exactly one"

    # Check for atoms other than carbon and oxygen
    allowed_atomic_nums = {6, 8}  # Carbon and Oxygen
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atomic_nums):
        return False, "Molecule contains atoms other than carbon and oxygen"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Molecule has {num_carbons} carbon atoms, fewer than required for a fatty acid"

    # Count number of oxygen atoms
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check that oxygen atoms are not excessive compared to carbon atoms
    if num_oxygens > num_carbons / 2:
        return False, f"Molecule has {num_oxygens} oxygen atoms, too many relative to carbon atoms"

    return True, "Molecule is a fatty acid anion with a carboxylate group and sufficient carbon chain"