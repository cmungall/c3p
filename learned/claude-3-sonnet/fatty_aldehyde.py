"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has an aldehyde group at one end of a carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde groups (-CH=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CH1,CH0]([H,#1,*])=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    if len(aldehyde_matches) > 1:
        return False, "Multiple aldehyde groups found"

    # Reject molecules with rings (except small rings that might be errors in SMILES)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        # Check ring sizes - only allow 3-membered rings which might be SMILES artifacts
        ssr = Chem.GetSymmSSSR(mol)
        for ring in ssr:
            if len(ring) > 3:
                return False, "Contains ring structure - not a fatty aldehyde"

    # Check for ester groups
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester group"

    # Check for carboxylic acids
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains carboxylic acid group"

    # Count carbons and check chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Allow smaller aldehydes if they have appropriate structure
    if carbon_count < 3:
        return False, "Carbon chain too short for fatty aldehyde"

    # Verify the aldehyde is terminal
    aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_matches[0][0])
    carbon_neighbors = [n for n in aldehyde_carbon.GetNeighbors() 
                       if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Aldehyde group is not terminal"

    # Check for excessive heteroatoms
    heteroatoms = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() not in (1, 6))
    if heteroatoms > carbon_count * 0.3:  # Allow some heteroatoms but not too many
        return False, "Too many heteroatoms for fatty aldehyde"

    # Check for connected carbon chain from aldehyde
    def trace_carbon_chain(atom, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                length = trace_carbon_chain(neighbor, visited)
                max_length = max(max_length, length)
        return max_length + 1

    chain_length = trace_carbon_chain(carbon_neighbors[0], set())
    if chain_length < 2:
        return False, "Insufficient carbon chain length"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    has_double_bonds = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bonds:
        return True, "Unsaturated fatty aldehyde with terminal aldehyde group"
    else:
        return True, "Saturated fatty aldehyde with terminal aldehyde group"