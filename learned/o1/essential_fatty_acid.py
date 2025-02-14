"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:18035 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid that cannot be synthesized
    by the human body and must be obtained from the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group (-C(=O)O)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Get the carbon chain attached to the carboxylic acid group
    # First, find the carboxylic acid carbon
    matches = mol.GetSubstructMatches(carboxylic_acid)
    carboxylic_carbons = [match[0] for match in matches]
    if not carboxylic_carbons:
        return False, "No carboxylic acid carbon found"

    carboxylic_carbon = carboxylic_carbons[0]
    
    # Perform a BFS to find the hydrocarbon chain
    visited = set()
    chain_atoms = []
    queue = [carboxylic_carbon]
    
    while queue:
        atom_idx = queue.pop(0)
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_atoms.append(atom_idx)
            # Add neighbors to queue
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() in (1, 6):  # Hydrogen or Carbon
                    queue.append(neighbor_idx)
                    
    # Count the number of carbons in the chain
    num_carbons = len(chain_atoms)
    if num_carbons < 10:
        return False, f"Carbon chain too short ({num_carbons} carbons), not a fatty acid"

    # Identify double bonds in the chain
    double_bond_positions = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Calculate position relative to carboxyl carbon
                try:
                    path_length = Chem.GetShortestPath(mol, carboxylic_carbon, begin_atom.GetIdx())
                    pos = len(path_length) - 1  # Position from carboxyl carbon
                    double_bond_positions.append(pos)
                except:
                    continue

    num_double_bonds = len(set(double_bond_positions))
    if num_double_bonds < 2:
        return False, f"Only {num_double_bonds} double bond(s) found, not polyunsaturated"

    # Check if any double bonds are beyond delta-9 position
    essential = any(pos > 9 for pos in double_bond_positions)
    if not essential:
        return False, "No double bonds beyond delta-9 position, not essential"

    return True, (
        f"Molecule is a fatty acid with {num_carbons} carbons, "
        f"{num_double_bonds} double bonds beyond delta-9 position"
    )