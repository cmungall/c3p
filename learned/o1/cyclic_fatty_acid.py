"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is any fatty acid containing anywhere in its structure a ring of atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    
    # Check for ring structures
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structures found"

    # Identify aliphatic chain connected to the carboxylic acid
    # We will look for a chain of at least 6 carbons attached to the carboxylic carbon
    fatty_acid_chain_length = 6  # Minimum chain length for fatty acids

    # Get indices of carboxylic acid carbon atoms
    carboxylic_acid_c_indices = [match[0] for match in carboxylic_acid_matches]

    # Check each carboxylic acid group for a long aliphatic chain
    for c_idx in carboxylic_acid_c_indices:
        # Use BFS to traverse the molecule from the carboxylic carbon
        visited = set()
        to_visit = [(c_idx, 0)]  # (atom index, current chain length)
        max_chain_length = 0
        
        while to_visit:
            current_idx, chain_length = to_visit.pop(0)
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            
            # Skip non-carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            
            # Update max chain length
            if chain_length > max_chain_length:
                max_chain_length = chain_length
            
            # Traverse neighbors
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(current_idx, n_idx)
                # Skip if already visited or if the bond is part of a ring
                if n_idx in visited or bond.IsInRing():
                    continue
                # Consider only single or double bonds (ignore triple bonds)
                if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
                    continue
                # Increment chain length if neighbor is carbon
                if neighbor.GetAtomicNum() == 6:
                    to_visit.append((n_idx, chain_length + 1))
        
        if max_chain_length >= fatty_acid_chain_length:
            return True, "Contains ring(s), carboxylic acid group, and long aliphatic chain"
    
    return False, f"No long aliphatic chain of at least {fatty_acid_chain_length} carbons found connected to carboxylic acid"