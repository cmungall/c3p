"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: straight-chain saturated fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Check for rings
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains cyclic structures"
        
    # Get the carboxylic acid carbon
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_matches[0][0])
    
    # Traverse the carbon chain from the carboxylic group
    visited = set()
    current = carboxyl_carbon
    chain_length = 0
    
    while current is not None:
        visited.add(current.GetIdx())
        chain_length += 1
        
        # Find next carbon in chain
        next_carbon = None
        for neighbor in current.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                # Check for branching
                carbon_neighbors = sum(1 for n in neighbor.GetNeighbors() 
                                    if n.GetAtomicNum() == 6)
                if carbon_neighbors > 2:
                    return False, "Branched carbon chain detected"
                next_carbon = neighbor
                break
        current = next_carbon

    # Check minimum chain length (typically 4+ carbons for fatty acids)
    if chain_length < 4:
        return False, "Carbon chain too short for a fatty acid"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains double bonds (unsaturated)"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains triple bonds (unsaturated)"

    # Check for non C,H,O atoms (excluding deuterium)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:
            if not (atomic_num == 1 and atom.GetIsotope() == 2):  # Allow deuterium
                return False, "Contains elements other than C, H, O (or deuterium)"

    # Check for ether linkages
    ether_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    if mol.HasSubstructMatch(ether_pattern):
        return False, "Contains ether linkages"

    # Check for additional oxygen-containing groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            # Skip the carboxylic acid oxygens
            if atom.GetIdx() not in [match[1] for match in carboxyl_matches] and \
               atom.GetIdx() not in [match[2] for match in carboxyl_matches]:
                return False, "Contains additional oxygen-containing groups"

    # All checks passed
    return True, "Straight-chain saturated fatty acid with single carboxylic acid group"