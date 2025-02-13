"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI:36942 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is any fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group(s)
    hydroxy_pattern = Chem.MolFromSmarts("O")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    # Look for aliphatic carbon chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No aliphatic carbon chain found"

    # Count rotatable bonds to verify aliphatic chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Carbon chain too short to be a fatty acid"

    # Check molecular weight - fatty acids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 6:
        return False, "Too few carbons for fatty acid"
    if o_count < 2:
        return False, "Not enough oxygens (must have at least carboxyl and hydroxy)"

    return True, "Contains a carboxylic acid group, aliphatic carbon chain, and at least one hydroxy group"