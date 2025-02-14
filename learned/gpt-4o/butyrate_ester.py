"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is defined as any carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define butyrate ester pattern specifically (butyric acid as ester: CCCC(=O)O)
    butyrate_esters_smarts = "CCCC(=O)O[*]"  # Finds esters where the esterifying group is bonded to the "O" atom
    butyrate_pattern = Chem.MolFromSmarts(butyrate_esters_smarts)
    if not mol.HasSubstructMatch(butyrate_pattern):
        return False, "No butyrate ester group found"

    return True, "Contains a butyrate ester group"