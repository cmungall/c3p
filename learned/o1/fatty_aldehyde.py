"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde (CHEBI:35581)
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check if molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), not a fatty aldehyde."
    
    # Check for disallowed elements (only C, H, O allowed)
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element with atomic number {atom.GetAtomicNum()}."
    
    # Define the aldehyde group pattern (terminal aldehyde)
    # Carbonyl carbon (C=O) connected to one carbon and implicit hydrogen
    aldehyde_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 2:
            has_double_bond_to_oxygen = False
            has_single_bond_to_carbon = False
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_bond_to_oxygen = True
                elif neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    has_single_bond_to_carbon = True
            if has_double_bond_to_oxygen and has_single_bond_to_carbon:
                aldehyde_carbons.append(atom)
    
    if not aldehyde_carbons:
        return False, "No terminal aldehyde group found."
    
    # Check if any aldehyde group is at the end of a sufficiently long carbon chain
    for aldehyde_c in aldehyde_carbons:
        # Get the carbon connected to the aldehyde carbon
        chain_atoms = set()
        for neighbor in aldehyde_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                chain_atom_idx = neighbor.GetIdx()
                break
        else:
            continue  # No carbon neighbor found, continue to next aldehyde carbon
        
        # Traverse the carbon chain starting from the chain_atom_idx
        visited = set()
        def traverse_chain(atom_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            visited.add(atom_idx)
            count = 1 if atom.GetAtomicNum() == 6 else 0
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                neighbor_atom = mol.GetAtomWithIdx(n_idx)
                if n_idx not in visited and neighbor_atom.GetAtomicNum() == 6:
                    count += traverse_chain(n_idx)
            return count
        
        chain_length = traverse_chain(chain_atom_idx)
        
        # Check if chain length is at least 4 carbons
        if chain_length >= 4:
            return True, "Molecule is a fatty aldehyde with a terminal aldehyde group and appropriate carbon chain."
        else:
            continue  # Try next aldehyde group
    
    return False, "No suitable terminal aldehyde group with sufficient carbon chain found."