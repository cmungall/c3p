"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene has at least one benzene ring with exactly three chlorine atoms attached to it
    directly.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)

    if not benzene_matches:
        return False, "No benzene ring found"

    # Iterate over each benzene ring match
    for match in benzene_matches:
        benzene_atoms = [mol.GetAtomWithIdx(i) for i in match]
        chlorine_count = 0
        for atom in benzene_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:
                    chlorine_count += 1
        if chlorine_count == 3:
           return True, "Contains a benzene ring with exactly 3 directly attached chlorine atoms"

    return False, "No benzene ring found with exactly 3 directly attached chlorine atoms"