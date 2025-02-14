"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    Requires both an L-alpha-amino acid backbone and an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for L-alpha-amino acid backbone:
    # R-[C@H](N)-C(=O)-O, where R can be any side chain
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)[C](=O)[O]")
    
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found"

    # Pattern for N-acyl substituent: N-C(=O)-R, where R is any group
    acyl_pattern = Chem.MolFromSmarts("N[C](=O)")

    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No N-acyl substituent found"
    
    return True, "Contains L-alpha-amino acid backbone with an N-acyl substituent"