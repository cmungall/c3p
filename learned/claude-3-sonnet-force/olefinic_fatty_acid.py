"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: CHEBI:37722 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is defined as any fatty acid containing at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for at least one C=C double bond
    olefin_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(olefin_pattern):
        return False, "No C=C double bonds found"
    
    # Check for long carbon chains (fatty acids are typically > 8 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Carbon chain too short for a fatty acid"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12:
        return False, "Too few carbons for a fatty acid"
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (carboxylic acid group)"
    
    return True, "Contains a carboxylic acid group and at least one C=C double bond"