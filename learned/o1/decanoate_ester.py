"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation
    of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester SMARTS pattern (carbonyl carbon connected to an ester oxygen)
    ester_smarts = '[C;X3](=O)[O;X2][#6]'
    ester_mol = Chem.MolFromSmarts(ester_smarts)
    esters = mol.GetSubstructMatches(ester_mol)

    if not esters:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain is decanoyl (10 carbons, linear, unbranched)
    for ester in esters:
        carbonyl_idx = ester[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)

        # Function to get acyl chain length
        def get_acyl_chain_length(carbonyl_atom):
            visited_atoms = set()
            chain_atoms = []

            def traverse(atom, prev_atom):
                # Check for cycles
                if atom.GetIdx() in visited_atoms:
                    return False  # Cycle detected

                visited_atoms.add(atom.GetIdx())

                if atom.GetAtomicNum() != 6:
                    return False  # Not a carbon atom

                # Get neighbors excluding hydrogens and previously visited atom
                neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1 and nbr.GetIdx() != prev_atom.GetIdx()]

                # Check for branching (more than one neighbor)
                if len(neighbors) > 1:
                    return False  # Branching detected

                chain_atoms.append(atom)

                if len(neighbors) == 0:
                    return True  # Reached end of chain

                # Proceed to next atom
                return traverse(neighbors[0], atom)

            # Start from the alpha carbon (carbon adjacent to carbonyl carbon)
            neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]

            if len(neighbors) != 1:
                return 0, False  # Should have exactly one carbon neighbor

            alpha_carbon = neighbors[0]

            if not traverse(alpha_carbon, carbonyl_atom):
                return 0, False  # Non-linear chain or cycle detected

            # Include the carbonyl carbon in the count
            total_carbons = len(chain_atoms) + 1

            return total_carbons, True

        chain_length, is_linear = get_acyl_chain_length(carbonyl_atom)
        if is_linear and chain_length == 10:
            return True, "Contains decanoate ester group"

    return False, "No decanoate ester groups found"