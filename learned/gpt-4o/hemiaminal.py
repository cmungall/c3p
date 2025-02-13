"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for hemiaminal (C bonded to OH and NH2)
    hemiaminal_pattern = Chem.MolFromSmarts("[CX4]([OH])([NH2])")
    
    # Check if the molecule matches the hemiaminal pattern
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains a hemiaminal motif: C with OH and NH2 attached"
    else:
        return False, "No hemiaminal motif found"

# Example usage:
# result, reason = is_hemiaminal("OC(N)CC")
# print(f"Is hemiaminal: {result}, Reason: {reason}")