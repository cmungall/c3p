"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a carbon bonded to both a hydroxyl (-OH) and an amino (-NH2 or substituted amine) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a hemiaminal: a carbon with single bond to OH and single bond to NX3
    hemiaminal_pattern = Chem.MolFromSmarts("[CX4]([OH1])[NX3]") #X3 is any atom with 3 bonds

    #Check if the pattern is present
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Found a carbon bonded to both a hydroxyl and an amine group via single bonds"
    else:
        return False, "No hemiaminal structure found"