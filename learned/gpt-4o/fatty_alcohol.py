"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary and secondary hydroxyl groups without aromatic bonds
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OH]")  # Hydroxyl on sp3 carbon
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No suitable sp3 carbon-bound hydroxyl group found"

    # Check that molecule is not primarily aromatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic structure"

    # Check for existence of long carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4;H2,H3][CX4;H2][CX4;H2]")  # At least three-carbon aliphatic chain
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No suitable aliphatic carbon chain found"

    # Count the number of carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Fatty alcohols need at least 3 carbons and 1 hydroxyl group
    if c_count < 3 or o_count < 1:
        return False, "Too few carbon atoms or hydroxyl groups for a fatty alcohol"

    # Ensure it is majorly aliphatic (considering branching and lack of complex rings)
    total_heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if c_count / total_heavy_atoms < 0.85:
        return False, f"Structure is not primarily aliphatic with {c_count} out of {total_heavy_atoms} heavy atoms being carbon"

    return True, f"Contains suitable alcohol group with {c_count} carbon atoms forming an aliphatic chain"