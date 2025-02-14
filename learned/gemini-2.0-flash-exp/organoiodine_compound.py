"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one iodine atom
    if not any(atom.GetAtomicNum() == 53 for atom in mol.GetAtoms()):
         return False, "Does not contain iodine"

    # Define the carbon-iodine bond pattern
    # ~ means any bond type, * is a wildcard for any atom.
    carbon_iodine_pattern = Chem.MolFromSmarts("[I]~*") 

    # Check if the molecule has a carbon-iodine bond
    if mol.HasSubstructMatch(carbon_iodine_pattern):
      #now ensure that the iodine is attached to a carbon, by iterating over the iodine atoms
        for atom in mol.GetAtoms():
          if atom.GetAtomicNum() == 53:
             for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                  return True, "Contains at least one carbon-iodine bond"
        return False, "Iodine not bonded to carbon"    
    else:
        return False, "Does not contain a carbon-iodine bond"