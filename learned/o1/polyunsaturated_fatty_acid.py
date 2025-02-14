"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is an unbranched aliphatic carboxylic acid containing more than one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify terminal carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    carb_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carb_matches:
        return False, "No terminal carboxylic acid group found"

    # Assume the first match is the terminal carboxylic acid
    carb_atom_idx = carb_matches[0][0]  # Index of the carbonyl carbon

    # Find alpha carbon (carbon connected to the carboxyl carbon)
    carboxyl_carbon = mol.GetAtomWithIdx(carb_atom_idx)
    alpha_carbons = [neighbor for neighbor in carboxyl_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6]
    if not alpha_carbons:
        return False, "No alpha carbon connected to carboxyl group"
    alpha_carbon_idx = alpha_carbons[0].GetIdx()

    # Recursive function to traverse the unbranched chain
    def traverse_chain(atom_idx, visited):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        chain_length = 1  # Count current atom
        double_bond_count = 0
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            n_idx = neighbor.GetIdx()
            if n_idx == carb_atom_idx:
                continue  # Skip back to carboxyl carbon
            # Avoid cycles and ensure unbranched chain
            if neighbor.GetAtomicNum() == 6 and n_idx not in visited:
                # Ensure neighbor carbon is not in a ring
                if not neighbor.IsInRing():
                    # Ensure the neighbor carbon is only connected to max 2 carbons (unbranched)
                    carbon_neighbors = [n for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 6]
                    if len(carbon_neighbors) <= 2:
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            bond_double = 1
                        else:
                            bond_double = 0
                        sub_chain_length, sub_double_bond_count = traverse_chain(n_idx, visited)
                        chain_length += sub_chain_length
                        double_bond_count += bond_double + sub_double_bond_count
                    else:
                        return chain_length, double_bond_count  # Branching detected, stop traversal
        return chain_length, double_bond_count

    # Start traversal from alpha carbon
    visited_atoms = set([carb_atom_idx])  # Exclude carboxyl carbon
    chain_length, double_bond_count = traverse_chain(alpha_carbon_idx, visited_atoms)

    # Check if there are more than one double bond
    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} double bonds in unbranched aliphatic chain"
    else:
        return False, f"Contains {double_bond_count} double bond(s); requires more than one"