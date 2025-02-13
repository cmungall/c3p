"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:38839 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds.

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
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count double bonds
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if n_double_bonds != 3:
        return False, f"Found {n_double_bonds} double bonds, need exactly 3"
    
    # Check for long carbon chain (>= 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for polyunsaturated (at least 3 double bonds)
    bond_iter = mol.GetBonds()
    unsaturation_count = sum(1 for bond in bond_iter if bond.GetBondType() == Chem.BondType.DOUBLE)
    if unsaturation_count < 3:
        return False, "Not polyunsaturated (need at least 3 double bonds)"
    
    return True, "Contains a carboxylic acid group and 3 double bonds in a long carbon chain"