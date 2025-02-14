"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid (SCFA) based on its SMILES string.
    An SCFA is defined as an aliphatic monocarboxylic acid with a chain length of 6 or less carbons.
    If any non-hydrocarbon substituent is present, the compound is not regarded as an SCFA.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Ensure molecule is aliphatic (no aromatic atoms)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms"
    
    # Ensure all atoms are C, H, or O
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_nums:
            return False, f"Contains atom other than C, H, O: {atom.GetSymbol()}"
    
    # Get atoms involved in the carboxylic acid group
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)
    
    # Count total oxygen atoms
    total_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if total_oxygen_atoms != 2:
        return False, f"Contains {total_oxygen_atoms} oxygen atoms, should have exactly 2 (from carboxylic acid group)"
    
    # Check for additional oxygen atoms not in carboxylic acid group
    oxygen_atoms_in_carboxylic = [idx for idx in carboxylic_acid_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
    if len(oxygen_atoms_in_carboxylic) != 2:
        return False, "Carboxylic acid group does not contain exactly 2 oxygen atoms"
    
    # Ensure there are no extra oxygen atoms
    if total_oxygen_atoms != len(oxygen_atoms_in_carboxylic):
        return False, "Contains additional oxygen atoms outside of carboxylic acid group"
    
    # Check for other disallowed functional groups (e.g., nitrogen, halogens, sulfur, phosphorus)
    disallowed_atomic_nums = {7, 9, 15, 16, 17, 35, 53}  # N, F, P, S, Cl, Br, I
    if any(atom.GetAtomicNum() in disallowed_atomic_nums for atom in mol.GetAtoms()):
        return False, "Contains disallowed heteroatoms"
    
    # Validate the carbon chain length (longest chain including the carboxylic carbon)
    # Identify the carboxylic acid carbon index
    carboxylic_carbon_idx = None
    for match in carboxylic_acid_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                carboxylic_carbon_idx = idx
                break
        if carboxylic_carbon_idx is not None:
            break
    if carboxylic_carbon_idx is None:
        return False, "Could not find carboxylic acid carbon atom"
    
    # Perform a depth-first search to find the longest carbon chain starting from carboxylic carbon
    visited = set()
    max_chain_length = [0]
    
    def dfs(node_idx, length):
        visited.add(node_idx)
        atom = mol.GetAtomWithIdx(node_idx)
        if atom.GetAtomicNum() != 6:
            # Only consider carbon atoms in the chain length
            pass
        else:
            if length > max_chain_length[0]:
                max_chain_length[0] = length
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(node_idx, n_idx)
                if n_idx not in visited and bond is not None:
                    n_atom = mol.GetAtomWithIdx(n_idx)
                    if n_atom.GetAtomicNum() == 6:
                        dfs(n_idx, length + 1)
        visited.remove(node_idx)
    
    dfs(carboxylic_carbon_idx, 1)
    
    if max_chain_length[0] > 6:
        return False, f"Carbon chain length is {max_chain_length[0]}, which is more than 6"
    
    return True, "Molecule is an aliphatic monocarboxylic acid with chain length of 6 or less and no additional functional groups"