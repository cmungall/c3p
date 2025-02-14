"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the 6-aminopurine moiety
    # 'c1nc2c(n1)ncn2[NH2]'
    # This specifies a 6-membered ring with 2 nitrogens, fused to a 5-membered ring with 2 nitrogens
    # and an amino group at position 6 of the 6 membered ring.
    
    adenine_smarts = Chem.MolFromSmarts("c1nc2c(n1)ncn2[NH2]")

    # Check for the presence of the 6-aminopurine substructure
    if mol.HasSubstructMatch(adenine_smarts):
        return True, "Contains a 6-aminopurine moiety"
    else:
        return False, "Does not contain a 6-aminopurine moiety"