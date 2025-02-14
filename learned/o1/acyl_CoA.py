"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:37577 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the Coenzyme A moiety (simplified for key features)
    coenzymeA_smarts = Chem.MolFromSmarts("""
    O[C@H]1[C@H](O)[C@@H](OP(=O)(O)O)[C@H](O[C@@H]1COP(=O)(O)OP(=O)(O)O)N2C=NC3=C2N=CN=C3N
    """)
    if not mol.HasSubstructMatch(coenzymeA_smarts):
        return False, "Coenzyme A moiety not found"

    # Define SMARTS pattern for the thioester linkage (S-C(=O)-C)
    thioester_smarts = Chem.MolFromSmarts("SC(=O)C")
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage not found"

    # Check for the acyl group attached via thioester bond
    # Extract the sulfur atom involved in thioester bond
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    else:
        # Optionally, verify that the acyl chain is attached to the sulfur via the thioester bond
        acyl_chain_length = 0
        sulfur_idx = thioester_matches[0][0]
        # Perform a BFS to count carbons attached to the carbonyl carbon (excluding CoA part)
        visited = set()
        queue = [thioester_matches[0][2]]  # Index of the carbon next to C=O (acyl chain)
        while queue:
            atom_idx = queue.pop(0)
            if atom_idx not in visited:
                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    acyl_chain_length += 1
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                            queue.append(neighbor_idx)

        if acyl_chain_length == 0:
            return False, "No acyl chain detected in thioester linkage"

    return True, "Contains Coenzyme A moiety linked via thioester bond to an acyl group"