"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered aromatic ring with three nitrogens at positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the 1,2,4-triazine core using aromatic nitrogens
    triazine_core_smarts = "n1cnncc1"
    triazine_pattern = Chem.MolFromSmarts(triazine_core_smarts)
    
    # Check if the molecule contains the core 1,2,4-triazine substructure
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "Molecule does not contain the 1,2,4-triazine core."

    return True, "Molecule contains the 1,2,4-triazine core."