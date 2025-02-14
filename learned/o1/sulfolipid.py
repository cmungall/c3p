"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:35346 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a 
    carbon-sulfur or oxygen-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted definition based on examples: sulfonic acid attached via O-S bond
    # Define SMARTS pattern for sulfonic acid residue connected via O-S bond
    sulfo_pattern = Chem.MolFromSmarts("[O;!$(*(=O))]-S(=O)(=O)-[O;H1,$([O-])]")  # Matches -O-S(=O)(=O)-O(H)
    if sulfo_pattern is None:
        return False, "Invalid SMARTS pattern for sulfonic acid residue"

    # Check for the sulfonic acid group connected via O-S bond
    sulfo_matches = mol.GetSubstructMatches(sulfo_pattern)
    if not sulfo_matches:
        return False, "No sulfonic acid residue connected via oxygen-sulfur bond found"

    # Detect long aliphatic chains (lipid chains)
    # Define a function to find aliphatic chains of length >=10
    def has_long_aliphatic_chain(mol, min_length=10):
        # Exclude rings and heteroatoms
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
        
        # Build adjacency list for carbons
        adj_list = {}
        for atom in carbon_atoms:
            idx = atom.GetIdx()
            adj_list[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                             if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
        
        # Use DFS to find chains
        visited = set()
        def dfs(current_atom, length):
            visited.add(current_atom)
            max_length = length
            for neighbor in adj_list.get(current_atom, []):
                if neighbor not in visited:
                    new_length = dfs(neighbor, length + 1)
                    if new_length > max_length:
                        max_length = new_length
            visited.remove(current_atom)
            return max_length

        # Check for chains starting from each carbon atom
        for atom in carbon_atoms:
            max_chain_length = dfs(atom.GetIdx(), 1)
            if max_chain_length >= min_length:
                return True
        return False

    # Check if the molecule has at least one long aliphatic chain
    has_lipid_chain = has_long_aliphatic_chain(mol, min_length=10)
    if not has_lipid_chain:
        return False, "No long aliphatic carbon chains (lipid) found"

    return True, "Contains sulfonic acid residue connected via oxygen-sulfur bond to a lipid"