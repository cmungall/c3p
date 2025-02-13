"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    'sphingomyelin d18:1' is defined as having sphingosine as the sphingoid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sphingosine backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("O[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@@H](/C=C/CCCCCCCCCCCCC)NC=O")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (d18:1) pattern found"

    # Define the amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    amide_count = len(mol.GetSubstructMatches(amide_pattern))
    if amide_count < 1:
        return False, f"Found {amide_count} amide linkages, but need at least 1"

    # Define the phosphocholine moiety pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine moiety found"

    return True, "Molecule has sphingosine backbone and phosphocholine group"