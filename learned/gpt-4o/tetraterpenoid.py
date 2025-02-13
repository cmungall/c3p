"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (35 <= c_count <= 55):
        return False, f"Expected around 40 carbon atoms, found {c_count}"

    # Check for an extended conjugated system typical in tetraterpenoids
    polyene_patterns = [
        Chem.MolFromSmarts("C=CC=CC=CC=CC=C"),  # Standard long conjugated chain
        Chem.MolFromSmarts("C=CC=C(C=C)C=C"),  # Branched conjugation
        Chem.MolFromSmarts("C=CC=CC=C(CC)C=C"),  # With terminal groups
    ]
    found_polyene = any(mol.HasSubstructMatch(pat) for pat in polyene_patterns)
    if not found_polyene:
        return False, "No extended conjugated polyene system found"

    # Check for presence of possible functional groups
    functional_patterns = [
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl group
        Chem.MolFromSmarts("[C]=O"),  # Carbonyl group
        Chem.MolFromSmarts("O([CX3]=[OX1])"),  # Ester or similar
        Chem.MolFromSmarts("O"),  # Any oxygen atom
    ]
    found_functional_groups = any(mol.HasSubstructMatch(pat) for pat in functional_patterns)

    if not found_functional_groups:
        return False, "Missing expected functional groups indicative of tetraterpenoids"
    
    return True, "Structure consistent with a tetraterpenoid: extended polyene and presence of functional groups"