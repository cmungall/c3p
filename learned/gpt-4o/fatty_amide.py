"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid,
    characterized by amide groups and usually long aliphatic hydrocarbon chains.

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

    # Look for amide group (primary, secondary, or tertiary)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check for long carbon chain or set of chains totaling a sufficient carbon count
    # Allow detection of chains with varied saturation or branching patterns
    long_chain_pattern1 = Chem.MolFromSmarts("C(CCCCCCCCC)")
    long_chain_pattern2 = Chem.MolFromSmarts("C=C(CCCCCC)")
    total_carbons = sum(atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() for atom in mol.GetAtoms())
    
    if len(mol.GetSubstructMatches(long_chain_pattern1)) > 0 or len(mol.GetSubstructMatches(long_chain_pattern2)) > 0:
        if total_carbons >= 12:
            return True, "Contains amide group and a sufficiently long aliphatic chain characteristic of fatty acids"

    return False, "Does not match typical characteristics of fatty amides or chain is too short"