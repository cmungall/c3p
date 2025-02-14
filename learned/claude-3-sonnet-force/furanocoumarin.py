"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:27727 furanocoumarin

A furanocoumarin is defined as any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coumarin substructure
    coumarin_pattern = Chem.MolFromSmarts("c1cc2ccoc2cc1=O")
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)

    # Look for furan substructure
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    furan_matches = mol.GetSubstructMatches(furan_pattern)

    # Check if both coumarin and furan rings are present and fused
    if coumarin_matches and furan_matches:
        # Check for atom sharing between coumarin and furan rings
        atom_idx_set = set()
        for match in coumarin_matches:
            atom_idx_set.update(match)
        for match in furan_matches:
            shared_atoms = atom_idx_set.intersection(match)
            if shared_atoms:
                return True, "Contains a fused coumarin and furan ring system"

    return False, "Does not contain a fused coumarin and furan ring system"