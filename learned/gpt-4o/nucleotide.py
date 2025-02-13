"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for nucleobases, sugar, and phosphate groups
    nucleobase_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2 | N-containing heterocycle")
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H]1 | Ribose or deoxyribose structure")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")  # Phosphate group

    # Check for nucleobase presence
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Check for sugar structure presence
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar structure (ribose or deoxyribose) found"
    
    # Check for phosphate group presence
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    return True, "Contains nucleobase attached to a sugar with a phosphate group"