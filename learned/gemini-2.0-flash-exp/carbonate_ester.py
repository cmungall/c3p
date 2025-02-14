"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the general structure R-O-C(=O)-O-R', where R and R' are organyl groups, 
    or R-O-C(=O)-O-H for a hydrogen carbonate. This also covers cyclic carbonates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbonate ester SMARTS pattern (non-cyclic, two organic groups):  [!H][O][CX3](=[OX1])[O][!H]
    carbonate_pattern = Chem.MolFromSmarts("[!H][O][C](=[O])([O][!H])")

    # Define SMARTS for cyclic carbonates, such as ethylene carbonate.
    cyclic_carbonate_pattern = Chem.MolFromSmarts("[O]1-[C](=[O])-[O]-[A]1")
    
    # Define SMARTS for hydrogen carbonates (R-O-C(=O)-OH)
    hydrogen_carbonate_pattern = Chem.MolFromSmarts("[!H][O][C](=[O])[OH1]")


    # Check if the molecule has the carbonate group using any of the patterns.
    if mol.HasSubstructMatch(carbonate_pattern):
        return True, "Contains a carbonate group with two organic groups"
    elif mol.HasSubstructMatch(cyclic_carbonate_pattern):
         return True, "Contains a cyclic carbonate group"
    elif mol.HasSubstructMatch(hydrogen_carbonate_pattern):
        return True, "Contains a hydrogen carbonate group"
    else:
      return False, "No carbonate group found"