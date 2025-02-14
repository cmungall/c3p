"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone and one acyl group (ester bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone: CH2(OH)-CH(OH)-CH2(OH)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for exactly one ester linkage: -O-C(=O)-
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 1"

    # Ensure the acyl group is sufficiently long (e.g., 8+ carbons in chain)
    # We will count carbon atoms starting from the carbonyl carbon in the ester to a terminal methyl
    acyl_carbon = ester_matches[0][1]  # The carbonyl carbon in the ester match
    carbon_count = 0
    # Use a set to track visited atoms to prevent cycles
    visited_atoms = set()
    
    # DFS or BFS could be used, but for simplicity and typical linear nature of fatty acids, iterate directly
    for atom in mol.GetAtomWithIdx(acyl_carbon).GetNeighbors():
        # Start counting from the carbon attached to the ester group
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in visited_atoms:
            # Move linearly along the chain
            current_atom = atom
            while current_atom.GetDegree() == 2 or current_atom.GetDegree() == 3:
                carbon_count += 1
                visited_atoms.add(current_atom.GetIdx())
                next_neighbors = [neighbor for neighbor in current_atom.GetNeighbors()
                                  if neighbor.GetIdx() != acyl_carbon and neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms]
                if next_neighbors:
                    current_atom = next_neighbors[0]
                else:
                    break

    # Check if the acyl chain length is >=8 carbons
    if carbon_count < 8:
        return False, f"Acyl chain length is {carbon_count}, which is too short for typical monoacylglycerol"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"