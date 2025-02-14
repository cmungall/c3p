"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxyl group on the carbon
    in the beta- or 3-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a hydroxy group at the beta position, next to carboxyl
    beta_hydroxy_pattern = Chem.MolFromSmarts("C(=O)C[C](O)")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group found at beta position"

    # Check for fatty acid chain length (at least 4 carbons, with one end = carboxylic acid) and at least 4 rotatable bonds to exclude small molecules.
    # We already checked that there is a carboxyl group, so we can check for at least 3 additional carbons for a total of 4
    chain_pattern = Chem.MolFromSmarts("C-C-C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Too short hydrocarbon chain"
    
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
         return False, "Too few rotatable bonds for fatty acid"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 4 or o_count < 3: # at least 4 carbons and 2 oxygens + the hydroxyl, for a total of 3
        return False, "Too few carbons or oxygens for fatty acid"

    return True, "3-Hydroxy fatty acid identified"