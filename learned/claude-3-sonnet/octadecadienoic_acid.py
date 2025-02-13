"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def get_main_chain(mol, carboxyl_carbon):
    """Helper function to find the main carbon chain starting from carboxyl group"""
    visited = set()
    best_chain = []
    
    def dfs(atom, current_chain):
        nonlocal best_chain
        visited.add(atom.GetIdx())
        current_chain.append(atom.GetIdx())
        
        # If this chain is longer than our best chain, update best chain
        if len(current_chain) > len(best_chain):
            best_chain = current_chain[:]
            
        # Explore neighbors
        for bond in atom.GetBonds():
            next_atom = bond.GetOtherAtom(atom)
            if next_atom.GetIdx() not in visited and next_atom.GetAtomicNum() == 6:
                dfs(next_atom, current_chain)
                
        visited.remove(atom.GetIdx())
        current_chain.pop()
    
    start_atom = mol.GetAtomWithIdx(carboxyl_carbon)
    dfs(start_atom, [])
    return best_chain

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for allowed atoms (C, H, O only)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {1, 6, 8}:
            return False, "Contains atoms other than C, H, and O"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    carboxyl_carbon = carboxyl_matches[0][0]
    
    # Get main chain
    main_chain = get_main_chain(mol, carboxyl_carbon)
    
    # Verify chain length
    if len(main_chain) != 18:
        return False, f"Main chain has {len(main_chain)} carbons, must be exactly 18"
    
    # Create a submolecule containing only the main chain
    chain_atoms = set(main_chain)
    chain_bonds = set()
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in chain_atoms and bond.GetEndAtomIdx() in chain_atoms:
            chain_bonds.add(bond.GetIdx())
    
    # Count double bonds in main chain
    double_bond_count = 0
    for bond_idx in chain_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_count += 1
    
    if double_bond_count != 2:
        return False, f"Contains {double_bond_count} C=C double bonds in main chain, must be exactly 2"
    
    # Check for branching in main chain
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        non_h_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
        chain_neighbors = [n for n in non_h_neighbors if n.GetIdx() in chain_atoms]
        
        # For carbons in main chain (except ends), should have exactly 2 carbon neighbors in chain
        if atom_idx != main_chain[0] and atom_idx != main_chain[-1]:
            carbon_chain_neighbors = [n for n in chain_neighbors if n.GetAtomicNum() == 6]
            if len(carbon_chain_neighbors) != 2:
                return False, "Not a straight chain structure"
    
    # Check for rings
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, must be acyclic"
    
    return True, "C18 straight-chain fatty acid with 2 C=C double bonds"