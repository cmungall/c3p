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
    
    # Check for phytosphingosine backbone (C18 chain with multiple hydroxyl groups and an amine)
    # SMILES pattern: A long-chain with three hydroxyls (-OH) and an amino group (-NH2)
    phytosphingosine_pattern = Chem.MolFromSmarts("[C@H]([C@H](O)[C@H](O)CCCCCCCCCCCCCC)[C@H](O)CO")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Check for N-acyl linkage: presence of C(=O)-N
    nacyl_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(nacyl_pattern):
        return False, "No N-acyl linkage found"

    return True, "Contains a phytosphingosine backbone with an N-acyl linkage"