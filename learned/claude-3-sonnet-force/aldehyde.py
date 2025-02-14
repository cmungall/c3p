"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:26718 aldehyde
An aldehyde is a compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde pattern (C(=O)H)
    aldehyde_pattern = Chem.MolFromSmarts("C(=O)[H]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Check if the aldehyde group is part of a ring
    # Aldehydes in rings are not considered aldehydes in this definition
    atom_ids = mol.GetSubstructMatches(aldehyde_pattern)[0]
    aldehyde_atom = mol.GetAtomWithIdx(atom_ids[0])
    if aldehyde_atom.IsInRing():
        return False, "Aldehyde group is part of a ring"

    # Check for other carbonyl groups
    # If there are other carbonyl groups, it's not an aldehyde
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) > 1:
        return False, "Contains multiple carbonyl groups"

    # Check for other heteroatoms attached to the carbonyl carbon
    # If there are other heteroatoms, it's not an aldehyde
    allowed_atoms = [8, 6, 1]  # O, C, H
    for atom in aldehyde_atom.GetNeighbors():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Carbonyl carbon has other heteroatoms attached"

    return True, "Contains an aldehyde group (C(=O)H)"