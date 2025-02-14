"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate has an acetate group attached to a phenol.

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

    # Check for phenyl ring and acetate directly attached using SMARTS
    phenyl_acetate_pattern = Chem.MolFromSmarts("CC(=O)Oc1ccccc1")
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
         return True, "Contains a phenyl ring with an acetate group directly attached via a phenolic oxygen"
    
    return False, "No phenyl acetate structure found"