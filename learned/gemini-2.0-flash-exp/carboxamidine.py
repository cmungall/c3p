"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2, where R can be H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check if it is guanidine C(=N)(N)N
    guanidine_pattern = Chem.MolFromSmarts("[C](=[N])([N])[N]")
    if mol.HasSubstructMatch(guanidine_pattern):
        return False, "Molecule is a Guanidine"
    
    # More specific carboxamidine pattern [C](=[N])[N]
    core_pattern = Chem.MolFromSmarts("[C](=[N])[N]")
    matches = mol.GetSubstructMatches(core_pattern)

    if not matches:
        return False, "Core carboxamidine structure C(=N)-N not found."

    return True, "Molecule matches the core structure of a carboxamidine"