"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be,
    but not necessarily, alcoholic. Aliphatic diols are also called glycols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hydroxy group pattern (OH group)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxy_pattern is None:
        return False, "Error in hydroxy SMARTS pattern"

    # Find all hydroxy groups in the molecule
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    hydroxy_oxygen_idxs = [match[0] for match in hydroxy_matches]

    # Define carboxylic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_oxygen_idxs = [match[1] for match in carboxylic_acid_matches]

    # Exclude hydroxy groups that are part of carboxylic acids
    non_carboxylic_hydroxy_oxygen_idxs = set(hydroxy_oxygen_idxs) - set(carboxylic_oxygen_idxs)

    num_hydroxy_groups = len(non_carboxylic_hydroxy_oxygen_idxs)

    # Check if the molecule contains exactly two hydroxy groups
    if num_hydroxy_groups == 2:
        return True, "Molecule contains exactly two hydroxy groups"
    else:
        return False, f"Molecule contains {num_hydroxy_groups} hydroxy group(s), diol requires exactly two"

    # Additional reasoning for why classification was made
    # Include the positions of the hydroxy groups if desired