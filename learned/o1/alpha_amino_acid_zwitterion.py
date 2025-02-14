"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57844 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is an amino acid zwitterion obtained by transfer of a proton from 
    the carboxy group to the amino group of any alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for alpha-amino-acid zwitterion
    # Alpha carbon connected to [NH3+], [C(=O)[O-]], and any side chain [*]
    zwitterion_pattern = Chem.MolFromSmarts("[C;H1;X4]([N+;H3])([*])[C](=O)[O-]")
    if zwitterion_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Search for the pattern in the molecule
    if mol.HasSubstructMatch(zwitterion_pattern):
        return True, "Contains alpha-amino-acid zwitterion structure"
    else:
        return False, "Alpha-amino-acid zwitterion pattern not found"