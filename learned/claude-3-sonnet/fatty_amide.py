"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:34974 fatty amide
A monocarboxylic acid amide derived from a fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (-C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Look for long carbon chain (>= 8 carbons) attached to amide
    fatty_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    fatty_chain_matches = mol.GetSubstructMatches(fatty_chain_pattern)
    if not fatty_chain_matches:
        return False, "No fatty acid chain found (>= 8 carbons)"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be fatty amide"

    # Check molecular weight - fatty amides typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for fatty amide"

    # Count carbons and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 10:
        return False, "Too few carbons for fatty amide"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (amide group)"

    return True, "Contains amide group with a long carbon chain (fatty acid)"