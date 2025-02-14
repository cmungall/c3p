"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an alpha-amino acid core where the carboxyl group
    is replaced by an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for alpha carbon with amino and carboxyl/ester groups
    alpha_amino_acid_ester_pattern = Chem.MolFromSmarts("[CX4]([NX3])([CX3](=[OX1])[OX2])")
    
    #check that the pattern is matched
    if not mol.HasSubstructMatch(alpha_amino_acid_ester_pattern):
         return False, "Molecule does not contain an alpha-amino acid with an ester group"
    
    return True, "Molecule contains an alpha-amino acid ester."