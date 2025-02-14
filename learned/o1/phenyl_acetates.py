"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation 
    of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the phenyl acetate pattern
    # The pattern matches an acetyl group connected via an ester oxygen to an aromatic ring
    phenyl_acetate_pattern = Chem.MolFromSmarts('[CH3][C](=O)[O][c]')

    # Check for substructure match
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Contains phenyl acetate group"
    else:
        return False, "No phenyl acetate group found"