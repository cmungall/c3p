"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group on the alpha carbon next to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha-amino acid pattern: [N][C](C(=O)O)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("N[C](C(=O)O)")
    
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Contains alpha-amino acid structure"
    else:
        return False, "No alpha-amino acid pattern found"