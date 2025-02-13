"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of the carboxy group of decanoic acid
    with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ester linkage pattern with a 10-carbon chain
    ester_and_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(ester_and_chain_pattern):
        return False, "No decanoate ester linkage found"

    return True, "Contains decanoate ester linkage"