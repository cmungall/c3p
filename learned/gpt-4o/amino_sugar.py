"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is defined as any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Look for sugar-like ring structures typical of pyranoses (6-membered rings)
    sugar_pattern = Chem.MolFromSmarts("C1[C@H]([O])[C@@H]([O])[C@H]([O])[C@H](O)O1")  # Pyranose ring
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No typical sugar-like ring structure found"

    # Look for amino groups, possibly replacing -OH groups in these structures
    amino_pattern = Chem.MolFromSmarts("[CX4;R][NX3;H2,H1,H0]")  # Carbon attached to an NHx in a ring
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group replacing a hydroxy group found"

    return True, "Contains a sugar-like ring structure with one or more hydroxyl groups replaced by amino groups"