"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:35638 oxo fatty acid

An oxo fatty acid is any fatty acid containing at least one aldehydic or ketonic group
in addition to the carboxylic acid group.

Typical features:
- Long aliphatic chain
- Terminal carboxylic acid group (-C(=O)O)
- At least one aldehydic (-C(=O)H) or ketonic (-C(=O)-) group
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for oxo fatty acid
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("[C;$(C(=O)O)]" +  # Terminal carboxyl group
                                                 "[C;$(C(=O)[C,H])]" +  # Oxo group (aldehyde or ketone)
                                                 "[C;$(CCC)]")  # Long aliphatic chain
    
    # Check if the molecule matches the pattern
    matches = mol.GetSubstructMatches(oxo_fatty_acid_pattern)
    
    if not matches:
        return False, "Molecule does not match the oxo fatty acid pattern"
    
    # Optional: Check for a minimum carbon chain length
    # carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # if carbon_chain_length < 6:
    #     return False, f"Carbon chain length ({carbon_chain_length}) too short for an oxo fatty acid"
    
    return True, "Molecule matches the structure of an oxo fatty acid"