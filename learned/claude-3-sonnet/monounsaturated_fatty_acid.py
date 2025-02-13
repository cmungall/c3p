"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: CHEBI:36975 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has one double or triple bond in the fatty acid chain
    and singly bonded carbon atoms in the rest of the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for one double or triple bond
    bond_types = [bond.GetBondType() for bond in mol.GetBonds()]
    double_bonds = bond_types.count(Chem.BondType.DOUBLE)
    triple_bonds = bond_types.count(Chem.BondType.TRIPLE)
    if double_bonds + triple_bonds != 1:
        return False, "Found more than one double/triple bond"
    
    # Check for singly bonded carbon chain
    chain_pattern = Chem.MolFromSmarts("[C;H3]-[C;H2]-[C;H2]~[C;H2]~[C;H2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No singly bonded carbon chain found"
    
    # Count carbon and hydrogen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    # Check for reasonable carbon count (>=4)
    if c_count < 4:
        return False, "Too few carbon atoms for fatty acid"
    
    # Check for reasonable hydrogen count (>=6)
    if h_count < 6:
        return False, "Too few hydrogen atoms for fatty acid"
    
    return True, "Contains one double/triple bond and a singly bonded carbon chain with a carboxylic acid group"