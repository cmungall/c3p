"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA resulting from the condensation of the thiol group of coenzyme A (CoA) with the carboxyl group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for CoA moiety (adenine moiety as a representative part)
    coa_pattern = Chem.MolFromSmarts("NC1=NC=NC2=C1N=CN=C2N")  # Adenine ring in CoA
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"

    # Define SMARTS pattern for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume the acyl chain is attached to the carbonyl carbon of the thioester
    # For each thioester linkage, analyze the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        sulfur_idx = match[2]      # Sulfur atom index

        # Get attached atoms to carbonyl carbon excluding the O and S
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        neighbors = [atom.GetIdx() for atom in carbonyl_c.GetNeighbors() if atom.GetIdx() not in match]
        if not neighbors:
            continue  # No attached atoms; skip
        # Traverse the acyl chain starting from the carbon attached to the carbonyl carbon
        acyl_chain_atoms = set()
        visited_atoms = set()
        stack = neighbors

        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            acyl_chain_atoms.add(atom_idx)
            # Add neighbors to stack excluding those already visited or involved in the thioester bond
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms and neighbor.GetAtomicNum() == 6 and neighbor_idx not in match:
                    stack.append(neighbor_idx)

        # Analyze acyl chain for branching
        is_branched = False
        for atom_idx in acyl_chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Count the number of carbon neighbors in the acyl chain
            n_chain_carbons = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in acyl_chain_atoms)
            if n_chain_carbons > 2:
                is_branched = True
                break  # Found branching

        if is_branched:
            return True, "Contains CoA moiety with branched-chain fatty acyl chain attached via thioester linkage"

    return False, "No branched-chain fatty acyl-CoA found"