"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is a glycerol backbone with three fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for esterified glycerol backbone
    glycerol_ester_pattern = Chem.MolFromSmarts("""
        [
            CH2  # First carbon of glycerol
            -
            O    # Connected to oxygen
            -
            C(=O)  # Connected to carbonyl
        ]
        -
        [
            CH   # Second carbon of glycerol
            -
            O
            -
            C(=O)
        ]
        -
        [
            CH2  # Third carbon of glycerol
            -
            O
            -
            C(=O)
        ]
    """)
    if glycerol_ester_pattern is None:
        return False, "Invalid glycerol ester pattern"

    # Check for glycerol backbone esterified
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No esterified glycerol backbone found"

    # Find ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups attached to glycerol, need exactly 3"

    # Check for long carbon chains attached via ester bonds (fatty acid chains)
    fatty_acid_chains = []
    for match in ester_matches:
        ester_atom_index = match[2]  # Index of the oxygen atom in the ester linkage
        # Get the carbon atom attached to the carbonyl carbon (start of the fatty acid chain)
        fatty_acid_start = mol.GetAtomWithIdx(ester_atom_index).GetNeighbors()[0]
        c_chain_length = 0
        visited_atoms = set()
        atoms_to_visit = [fatty_acid_start.GetIdx()]

        # Traverse the fatty acid chain
        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                c_chain_length += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited_atoms and neighbor.GetAtomicNum() in [6, 1]:  # Carbon or Hydrogen
                        atoms_to_visit.append(neighbor_idx)

        if c_chain_length < 4:
            return False, "Fatty acid chains too short, need chains longer than 4 carbons"
        fatty_acid_chains.append(c_chain_length)

    # If all checks pass
    return True, "Contains esterified glycerol backbone with three fatty acid chains"