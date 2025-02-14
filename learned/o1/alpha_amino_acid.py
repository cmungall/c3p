"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is defined as an amino acid where the amino group is located
    on the carbon atom alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid
    # Amino group connected to alpha carbon, which is connected to a carboxyl group
    alpha_amino_acid_smarts = '[NX3;H2,H1;+0]-[CX4H]-[*]-[CX3](=O)[O;H1,-1]'

    pattern = Chem.MolFromSmarts(alpha_amino_acid_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the alpha-amino acid pattern
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Contains alpha-amino acid structure"
    else:
        return False, "No alpha-amino acid structure found"