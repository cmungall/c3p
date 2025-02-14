"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is defined as an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Count carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)

    # Count other atoms to determine the "fatty" nature
    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]]
    num_other = len(other_atoms)
    if num_other > 0:
      return False, "Contains atoms other than C, H, O"

    # Check if number of carbon is within the range
    if num_carbons < 3 :
        return False, f"Too few carbon atoms ({num_carbons}), must be at least 3"
    if num_carbons > 27 :
        #Check that it has a long carbon chain
         carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
         if not mol.HasSubstructMatch(carbon_chain_pattern):
              return False, "Has >27 carbons, but no long chain"

    
    if num_carbons >= 3 and num_carbons <=27:
         carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
         if not mol.HasSubstructMatch(carbon_chain_pattern):
              return False, "Must contain at least 3 carbons in an aliphatic chain"

    return True, "Meets criteria for a fatty alcohol (3 to >27 C, at least one OH)"