"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    Very long-chain fatty acyl-CoAs have fatty acid chains longer than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: (is_vlcfa_coa, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety
    # Look for adenine
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"

    # Look for phosphate groups (need at least 3 for CoA)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 3:
        return False, "Missing phosphate groups characteristic of CoA"

    # Count carbons in the main chain
    # First get the carbon atoms in the thioester group
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    thioester_carbon = thioester_matches[0][0]  # First carbon of the thioester group
    
    # Get all carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Find the longest carbon chain starting from the thioester carbon
    visited = set()
    def dfs_carbon_chain(atom_idx, depth=0):
        if atom_idx not in carbon_atoms or atom_idx in visited:
            return depth
        visited.add(atom_idx)
        max_depth = depth
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Only follow carbon atoms
                max_depth = max(max_depth, dfs_carbon_chain(neighbor.GetIdx(), depth + 1))
        return max_depth

    # Get chain length starting from thioester carbon
    chain_length = dfs_carbon_chain(thioester_carbon)
    
    if chain_length <= 22:
        return False, f"Fatty acid chain length (C{chain_length}) not greater than C22"
    
    if chain_length > 40:
        return False, f"Fatty acid chain length (C{chain_length}) unreasonably long"

    # Check for characteristic functional groups that may be present
    modifications = []
    
    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("CC=CC")
    if mol.HasSubstructMatch(double_bond_pattern):
        modifications.append("unsaturated")
    
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("CC(O)C")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        modifications.append("hydroxylated")
    
    # Check for oxo groups
    oxo_pattern = Chem.MolFromSmarts("CC(=O)C")
    if mol.HasSubstructMatch(oxo_pattern):
        modifications.append("oxo")

    modification_str = " and ".join(modifications)
    if modification_str:
        modification_str = f" ({modification_str})"
    
    return True, f"Very long-chain fatty acyl-CoA with C{chain_length} fatty acid chain{modification_str}"