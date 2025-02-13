"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:34842 organic sulfide
Organic sulfides are compounds with the structure RSR (R =/= H).
They were once called thioethers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check that each sulfur has two R groups attached (R =/= H)
    for s_atom in sulfur_atoms:
        r_groups = [nb.GetSymbol() != "H" for nb in s_atom.GetNeighbors()]
        if sum(r_groups) != 2:
            return False, f"Sulfur atom has {sum(r_groups)} non-hydrogen substituents, should be 2"

    # If all checks pass, it is an organic sulfide
    return True, "Contains sulfur atoms with two non-hydrogen substituents (RSR)"