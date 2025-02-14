"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:18010 nucleoside
An N-glycosyl compound that has both a nucleobase, normally adenine, guanine, xanthine, thymine, cytosine or uracil, and either a ribose or deoxyribose as functional parents.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ribose or deoxyribose sugar
    sugar_pattern = Chem.MolFromSmarts("[OX2]C[CX4]([OX2])[CX4]([CX4]([OX2])[CX4][OX2])[CH2][OX2]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Look for common nucleobases
    nucleobase_patterns = [
        Chem.MolFromSmarts("[$(nc1ncnc1N)]"),  # Adenine
        Chem.MolFromSmarts("[$(nc1nc(=O)nc(N)n1)]"),  # Guanine
        Chem.MolFromSmarts("[$(nc1nc(=O)c2ncn(n2)c1=O)]"),  # Xanthine
        Chem.MolFromSmarts("[$(Cc1c(N)nc(=O)[nH]c1=O)]"),  # Thymine
        Chem.MolFromSmarts("[$(NC(=O)nc1nc(N)ccn1)]"),  # Cytosine
        Chem.MolFromSmarts("[$(Oc1ccnc(=O)[nH]1)]")  # Uracil
    ]

    nucleobase_match = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not nucleobase_match:
        return False, "No common nucleobase found"

    # Check for N-glycosidic bond between sugar and nucleobase
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][NX3]")
    if not AllChem.FindMolChiralCenters(mol, glycosidic_bond_pattern):
        return False, "No N-glycosidic bond found between sugar and nucleobase"

    return True, "Contains a nucleobase (adenine, guanine, xanthine, thymine, cytosine or uracil) N-glycosidically linked to a ribose or deoxyribose sugar"