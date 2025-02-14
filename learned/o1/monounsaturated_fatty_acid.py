"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has a carboxylic acid group and a long carbon chain
    with exactly one double or triple bond, and singly bonded carbon atoms in the rest of the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Get the indices of the carboxylic acid carbon
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    carboxy_carbons = [match[0] for match in carboxy_matches]
    if not carboxy_carbons:
        return False, "No carboxylic acid carbon found"
    carboxy_carbon_idx = carboxy_carbons[0]

    # Identify the longest carbon chain starting from the carboxylic acid carbon
    def get_longest_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return []
        visited.add(atom_idx)
        max_chain = []
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                chain = get_longest_chain(neighbor_idx, visited.copy())
                if len(chain) > len(max_chain):
                    max_chain = chain
        return [atom_idx] + max_chain

    longest_chain = get_longest_chain(carboxy_carbon_idx, set())
    if len(longest_chain) < 2:
        return False, "Chain is too short to be a fatty acid"

    # Count double and triple bonds in the chain
    unsaturation_count = 0
    for i in range(len(longest_chain) - 1):
        bond = mol.GetBondBetweenAtoms(longest_chain[i], longest_chain[i+1])
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE or bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            unsaturation_count += 1
            bond_type = bond.GetBondType()
        elif bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, "Non-single bond found that is not double or triple"

    if unsaturation_count == 0:
        return False, "No double or triple bonds found in the chain"
    elif unsaturation_count > 1:
        return False, f"More than one double or triple bond found in the chain ({unsaturation_count} unsaturations)"

    # Check that other bonds are single bonds
    for i in range(len(longest_chain) - 1):
        bond = mol.GetBondBetweenAtoms(longest_chain[i], longest_chain[i+1])
        if bond.GetBondType() not in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            return False, "Bond type not single, double, or triple found in chain"

    # Ensure rest of the chain is saturated
    if unsaturation_count == 1:
        return True, "Molecule is a monounsaturated fatty acid"
    else:
        return False, "Unsaturation count does not match monounsaturated fatty acid definition"