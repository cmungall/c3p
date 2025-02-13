"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine consists of a phytosphingosine backbone with a fatty acyl group
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible pattern for the phytosphingosine backbone
    # This pattern covers the generic structure with three hydroxyl groups and an amine
    phytosphingosine_pattern = Chem.MolFromSmarts("O[C@H]([C@H](O)[C@H](CN)O)CCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Define a flexible pattern for an N-acyl linkage: presence of C(=O)-N
    nacyl_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(nacyl_pattern):
        return False, "No N-acyl linkage found"

    return True, "Contains a phytosphingosine backbone with an N-acyl linkage"