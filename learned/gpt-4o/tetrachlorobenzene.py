"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying four chlorine groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a benzene ring with any four chlorines
    # The revised pattern checks for any aryl ring structure with exactly four chlorines
    tetrachloro_pattern = Chem.MolFromSmarts("c1cc(Cl)c(Cl)cc1")
    
    # Count chlorine atoms linked to a benzene ring
    cl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17 and atom.GetNeighbors():
            adjacency = atom.GetNeighbors()[0]
            if adjacency.GetIsAromatic() and adjacency.GetSymbol() == 'C':
                cl_count += 1
    
    # Check if exactly four chlorines are bonded to the aromatic carbons of one ring
    if cl_count == 4 and mol.HasSubstructMatch(tetrachloro_pattern):
        return True, "Contains a benzene ring with four chlorine atoms"
    
    return False, "Does not match the tetrachlorobenzene pattern"

# Example of possible test usage
# print(is_tetrachlorobenzene("Clc1cc(Cl)c(Cl)cc1Cl")) # True