"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35519 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine contains the functional group RC(=NR)NR2, where a carbon atom
    is double-bonded to one nitrogen and single-bonded to another nitrogen.

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

    # Define the carboxamidine SMARTS pattern
    # Carbon (sp2) double-bonded to nitrogen and single-bonded to another nitrogen
    carboxamidine_pattern = Chem.MolFromSmarts('[CX3](=[NX2])[NX3;H1,H2,H0]')

    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine functional group"
    else:
        return False, "Does not contain carboxamidine functional group"