"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol where at least one hydroxyl group is attached to a non-aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one alcohol group (-OH) bound to an aliphatic carbon.
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    
    if not mol.HasSubstructMatch(alcohol_pattern):
      return False, "No aliphatic alcohol group found"

    #Check the matched carbon atoms to verify they are not aromatic
    for match in mol.GetSubstructMatches(alcohol_pattern):
        carbon_idx = match[0]
        c_atom = mol.GetAtomWithIdx(carbon_idx)
        if not c_atom.GetIsAromatic():
             return True, "Molecule contains at least one alcohol group attached to a non-aromatic carbon."
    
    return False, "No aliphatic alcohol group found"