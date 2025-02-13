"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thiol group pattern, -SH attached to carbon (aliphatic or aromatic)
    # This includes aliphatic carbon ([CX4H,CX3H0]), aromatic carbon ([c]), and covers general aliphatic connections as well.
    thiol_pattern_aliphatic = Chem.MolFromSmarts("[CX4H0,CX3][SH]")
    thiol_pattern_aromatic = Chem.MolFromSmarts("[c][SH]")

    # Check for thiol pattern matches in the molecule
    if mol.HasSubstructMatch(thiol_pattern_aliphatic) or mol.HasSubstructMatch(thiol_pattern_aromatic):
        return True, "Contains a thiol group (-SH) attached to a carbon atom"

    return False, "No thiol group found"