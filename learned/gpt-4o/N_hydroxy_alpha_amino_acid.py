"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    An N-hydroxy-alpha-amino-acid is defined as any amino acid in which at least one hydrogen 
    attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for N-hydroxy group specifically on an alpha carbon attached nitrogen
    n_hydroxy_alpha_amino_pattern = Chem.MolFromSmarts("[C;H1,H2](N(O))[C](=O)O")
    n_hydroxy_dihydroxy_alpha_amino_pattern = Chem.MolFromSmarts("[C;H1,H2](N(O)O)[C](=O)O")
    
    if mol.HasSubstructMatch(n_hydroxy_alpha_amino_pattern) or \
       mol.HasSubstructMatch(n_hydroxy_dihydroxy_alpha_amino_pattern):
        return True, "Contains N-hydroxy substitution on alpha-amino group and carboxylic acid characteristic of alpha-amino acids"

    return False, "Does not meet the criteria for N-hydroxy-alpha-amino-acid"