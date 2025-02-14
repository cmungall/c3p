"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a fatty acid with more than one carbon-carbon double or triple bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a long carbon chain attached to the acid group
    # This allows for various lengths and includes possible single bonds, double bonds, and single-bonded carbons.
    chain_pattern = Chem.MolFromSmarts("[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]~[CH2,CH1,CH0]")
    chain_match = mol.GetSubstructMatches(chain_pattern)
    if not chain_match:
        return False, "No long carbon chain with acid group found."

    # Check for more than one double or triple bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    triple_bond_matches = mol.GetSubstructMatches(triple_bond_pattern)
    if len(double_bond_matches) + len(triple_bond_matches) <= 1:
        return False, f"Found {len(double_bond_matches)} double and {len(triple_bond_matches)} triple bonds, need more than one."


    # Check that double bonds are part of the main fatty acid chain
    chain_double_bond = Chem.MolFromSmarts("[CH2,CH1,CH0]~[CH2,CH1,CH0]=[CH,CH0]~[CH2,CH1,CH0]")
    if not mol.HasSubstructMatch(chain_double_bond):
        return False, "Double/triple bonds not part of the main fatty acid chain."
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chain too short to be a fatty acid"

     # Count total carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons for a fatty acid"

    # Check molecular weight - PUFA typically between 200-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, "Molecular weight is outside the typical range of a PUFA"
        
    return True, "Contains a carboxylic acid group and more than one double/triple bond within the main fatty acid chain."