"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group 
    and a long hydrocarbon chain, possibly with some functionalization.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group pattern [CX3](=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found (required for fatty acid anion)"
    
    # Iterate over all carbon atoms to find longest carbon chain, considering
    # branches, rings, and functional groups:
    def explore_chain(atom, visited):
        if atom.GetAtomicNum() != 6 or atom.GetIdx() in visited:
            return 0
        visited.add(atom.GetIdx())
        max_chain_length = 1
        for neighbor in atom.GetNeighbors():
            max_chain_length = max(max_chain_length, 1 + explore_chain(neighbor, visited))
        visited.remove(atom.GetIdx())
        return max_chain_length

    longest_chain = max(explore_chain(atom, set()) for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if longest_chain < 12:  # At least 12 carbons in the longest chain
        return False, "Longest carbon chain is too short to be a fatty acid anion"

    # Optionally check the balance of functional groups vs. carbon atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Ratio for functional groups to carbon atoms:
    if o_count > c_count / 4:
        return False, "Too many functional groups for a typical fatty acid anion structure"
    
    return True, "Contains a carboxylate group and a long hydrocarbon chain characteristic of a fatty acid anion"