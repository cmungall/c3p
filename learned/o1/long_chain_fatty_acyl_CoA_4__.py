"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:77989 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Checks for the presence of a CoA moiety, a thioester linkage, and a long acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple:
            bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
            str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA SMARTS pattern (simplified)
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n2cnc3c(ncnc23)N"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    # Check for CoA moiety
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "CoA moiety not found"

    # Get atom indices of CoA moiety
    coa_atom_indices = set()
    for match in coa_matches:
        coa_atom_indices.update(match)

    # Define thioester linkage pattern: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume the first thioester match is the one connecting the acyl chain to CoA
    thioester_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon index

    # Traverse the acyl chain from the carbonyl carbon
    acyl_chain_length = 0
    visited_atoms = set()
    atoms_to_visit = [thioester_carbon_idx]

    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms or atom_idx in coa_atom_indices:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Increment chain length if not the carbonyl carbon
            if atom_idx != thioester_carbon_idx:
                acyl_chain_length += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms and neighbor_idx not in coa_atom_indices:
                    neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                    # Exclude sulfur atom (S) and limit traversal to carbon and hydrogen atoms
                    if neighbor_atom.GetAtomicNum() in [6, 1]:  # Carbon or Hydrogen
                        atoms_to_visit.append(neighbor_idx)

    if acyl_chain_length < 12:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons), need at least 12"

    return True, f"Contains CoA moiety with long acyl chain of {acyl_chain_length} carbons"