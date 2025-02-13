"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at the
    opposite end (omega position) of a primarily linear carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Look for terminal hydroxyl group (CH2-OH)
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2][OX2H1]")
    if not mol.HasSubstructMatch(terminal_oh_pattern):
        return False, "No terminal hydroxyl group found"

    # Get the indices of key atoms
    carboxylic_carbon = mol.GetSubstructMatch(carboxylic_pattern)[0]
    terminal_ch2 = mol.GetSubstructMatch(terminal_oh_pattern)[0]
    
    # Get shortest path between COOH carbon and terminal CH2
    path = Chem.GetShortestPath(mol, carboxylic_carbon, terminal_ch2)
    if not path:
        return False, "No valid path between COOH and terminal OH"
    
    # Check path length (minimum 4 atoms for shortest omega-hydroxy fatty acid)
    if len(path) < 4:
        return False, "Chain too short for omega-hydroxy fatty acid"
        
    # Check for excessive hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches > 3:  # Allow up to 3 OH groups (COOH + terminal OH + one possible extra)
        return False, "Too many hydroxyl groups for omega-hydroxy fatty acid"

    # Verify the atoms in the main chain
    path_atoms = [mol.GetAtomWithIdx(i) for i in path]
    non_carbon_count = sum(1 for atom in path_atoms if atom.GetAtomicNum() != 6)
    if non_carbon_count > 2:  # Allow only COOH and terminal OH as non-carbon
        return False, "Main chain contains too many non-carbon atoms"
    
    # Check for branching along the main chain
    for atom_idx in path:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Only check carbon atoms
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
                return False, "Excessive branching along main chain"

    # Count total carbons and check ratio against path length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    path_length = len(path)
    if carbon_count > path_length + 1:  # Allow at most one carbon branch
        return False, "Too many carbons outside main chain"

    # Verify the molecule consists only of C, H, O
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains elements other than C, H, O"

    return True, "Contains terminal hydroxyl and carboxylic acid groups at opposite ends of primarily linear carbon chain"