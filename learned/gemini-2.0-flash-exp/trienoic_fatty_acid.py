"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a fatty acid containing a carboxylic acid group and exactly three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly three double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 3:
      return False, f"Found {len(double_bond_matches)} double bonds, need exactly 3"
      
    # Check for long carbon chain (fatty acid part) - at least 8 carbons outside the carboxyl group
    
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_matches) < 1:
        return False, "No long carbon chain found"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 10 :
        return False, "Too few carbons for fatty acid"
    if o_count < 2:
      return False, "Not enough oxygens"

    # Check molecular weight - fatty acids typically > 150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for fatty acid"
    
    
    return True, "Contains carboxylic acid and exactly three double bonds and long carbon chain"