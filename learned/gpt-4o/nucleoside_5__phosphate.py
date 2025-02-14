"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribosyl or deoxyribosyl sugar linked to
    a nucleobase with a phosphate group attached to the 5' position of the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for ribose, deoxyribose, nucleobases, and phosphate groups
    ribose_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H](O)[C@H]1O)")
    deoxyribose_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H]([C@H]1)O)")
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(O)=O")

    # Check for the ribose or deoxyribose part
    if not (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
        return False, "No ribose or deoxyribose sugar detected"
    
    # Simplified check for presence of a common nucleobase (using canonical SMILES patterns)
    adenine_base = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    guanine_base = Chem.MolFromSmarts("c1[nH]c2c(n1)ncnc2")
    cytosine_base = Chem.MolFromSmarts("n1cnc(c=O)n1")
    thymine_base = Chem.MolFromSmarts("c1c[nH]c(=O)nc1=O")
    uracil_base = Chem.MolFromSmarts("c1c[nH]c(=O)n1")
    
    if not any(
        mol.HasSubstructMatch(base)
        for base in [adenine_base, guanine_base, cytosine_base, thymine_base, uracil_base]
    ):
        return False, "No recognizable nucleobase found"

    # Check for the phosphate group at 5' position
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group at the 5' position"

    return True, "Contains ribose or deoxyribose sugar, nucleobase, and 5'-phosphate group"