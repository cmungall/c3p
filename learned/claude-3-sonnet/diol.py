"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:27660 diol
A compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
Aliphatic diols are also called glycols.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound containing two hydroxy (-OH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))

    # Diol must have exactly 2 hydroxy groups
    if hydroxy_count != 2:
        return False, f"Found {hydroxy_count} hydroxy groups, need exactly 2"

    # Check for molecular weight range (typically < 500 Da for small diols)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for small diol"

    # Check for restricted functional groups
    # (diols should not contain other reactive groups like aldehydes, ketones, etc.)
    restricted_pattern = Chem.MolFromSmarts("[C$(C=O)][OX2H,OX1]")
    if mol.HasSubstructMatch(restricted_pattern):
        return False, "Contains restricted functional groups (aldehyde, ketone)"

    return True, "Contains exactly 2 hydroxy groups"