"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: carbonate ester (Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups.)
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has a central carbon atom doubly bonded to an oxygen atom
    and singly bonded to two oxygens, each attached to a carbon (ester groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a carbonate ester functional group
    # Pattern: carbon atom with double bond to oxygen and single bonds to two oxygens
    # Each oxygen is singly bonded to a carbon atom (i.e., O-C)
    carbonate_pattern = Chem.MolFromSmarts("[CX3](=O)(O[CX4,CX3,CX2])(O[CX4,CX3,CX2])")
    
    # Search for the carbonate ester pattern in the molecule
    if mol.HasSubstructMatch(carbonate_pattern):
        return True, "Contains carbonate ester functional group"
    else:
        return False, "Does not contain carbonate ester functional group"