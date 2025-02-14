"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:36976 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the longest chain length
    longest_chain = AllChem.InChI(mol).split('/')[1]
    chain_length = int(longest_chain.split('.')[0])
    
    # Check chain length
    if chain_length <= 22:
        return False, f"Chain length {chain_length} is too short for very long-chain fatty acid"
    
    # Check for carboxylic acid group (-C(=O)O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-,N]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for linear chain
    linear_pattern = Chem.MolFromSmarts("[CH3][CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(linear_pattern):
        return False, "Not a linear chain"

    # Check for additional functional groups
    if Chem.MolFromSmiles(smiles).GetRingInfo().NumRings() > 0:
        return False, "Contains rings, not a pure fatty acid"

    # Count carbons and hydrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if c_count < chain_length:
        return False, "Too few carbons for given chain length"
    if h_count < 2 * chain_length + 2:
        return False, "Too few hydrogens for linear chain"

    return True, f"Linear fatty acid with chain length {chain_length} (> C22)"