"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: octanoate ester
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the carboxylic acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2H0][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define octanoyl chain pattern (C(=O)CCCCCCC)
    octanoyl_pattern = Chem.MolFromSmarts("[C](=O)[CH2][CH2][CH2][CH2][CH2][CH2][CH3]")

    # Check each ester group for octanoyl chain
    for match in ester_matches:
        # Extract atoms from the match
        r1_c, c_carb, o_ester, r2_c = match[0], match[1], match[3], match[4]
        
        # Create a mapping for the pattern match
        mapping = {0: c_carb, 1: match[2], 2: match[3], 3: r2_c}
        
        # Get the sub-molecule representing the acyl chain
        acyl_chain = Chem.PathToSubmol(mol, [c_carb] + list(Chem.FindAtomEnvironmentOfRadiusN(mol, 7, c_carb)))
        
        # Check if the acyl chain matches the octanoyl pattern
        if acyl_chain.HasSubstructMatch(octanoyl_pattern):
            return True, "Contains octanoyl ester group"
        
        # Alternatively, traverse the acyl chain
        chain_atoms = set()
        stack = [(c_carb, None)]  # (current_atom, previous_atom)
        while stack:
            current_atom, previous_atom = stack.pop()
            atom = mol.GetAtomWithIdx(current_atom)
            if atom.GetAtomicNum() != 6:
                break  # Non-carbon atom found, not octanoyl
            chain_atoms.add(current_atom)
            neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous_atom]
            # Exclude oxygen atoms (carbonyl oxygen and ester oxygen)
            neighbors = [n for n in neighbors if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
            if len(chain_atoms) > 8:
                break  # Chain too long, not octanoyl
            if neighbors:
                stack.append((neighbors[0], current_atom))
            else:
                if len(chain_atoms) == 8:
                    return True, "Contains octanoyl ester group"
                break  # Reached end of chain

    return False, "No octanoyl ester groups found"