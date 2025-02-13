"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:36976 long-chain fatty acid
A fatty acid with a chain length ranging from C13 to C22.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Long-chain fatty acids have 13 to 22 carbon atoms
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Molecule has {carbon_count} carbon atoms, should be between 13 and 22"
    
    # Look for carboxylic acid group (-C(=O)OH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for long carbon chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "No long carbon chain found"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chain too short to be a long-chain fatty acid"
    
    return True, "Contains a carboxylic acid group and a long carbon chain (C13 to C22)"