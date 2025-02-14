"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:38165 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxy group at the 3rd carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chain (fatty acid)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "No long carbon chain found"
    
    # Look for hydroxyl group at 3rd carbon from carboxylic acid
    hydroxy_pattern = Chem.MolFromSmarts("[CX4]([CX3]([OX2H1])[CX3](=O)[OX2H0])[CX4,CX3]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, "No hydroxyl group at 3rd carbon found"
    
    return True, "Contains a carboxylic acid group and a hydroxyl group at the 3rd carbon"

# Example usage
smiles = "CCCCCCCCCCCCCCCCCCCCCCC[C@H](O)CC(O)=O"
is_3hydroxy, reason = is_3_hydroxy_fatty_acid(smiles)
print(is_3hydroxy, reason)