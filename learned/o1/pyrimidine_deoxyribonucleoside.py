"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base attached to a deoxyribose sugar via a β-N1-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group (to exclude nucleotides)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if phosphate_pattern is None:
        return None, "Invalid phosphate SMARTS pattern"
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Molecule is a nucleotide (contains phosphate group)"

    # Find pyrimidine ring
    pyrimidine_pattern = Chem.MolFromSmarts("n1cnccc1")
    if pyrimidine_pattern is None:
        return None, "Invalid pyrimidine SMARTS pattern"
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)
    if not pyrimidine_matches:
        return False, "No pyrimidine base found"

    # Find sugar ring (five-membered ring containing oxygen)
    sugar_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1")
    if sugar_pattern is None:
        return None, "Invalid sugar SMARTS pattern"
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No deoxyribose sugar found"

    # Now check for glycosidic bond between N1 of pyrimidine and C1' of sugar
    # For each pyrimidine ring
    for pyrimidine in pyrimidine_matches:
        n1_idx = pyrimidine[0]  # N1 atom index
        n1_atom = mol.GetAtomWithIdx(n1_idx)
        # Check bonds from N1
        for bond in n1_atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(n1_atom)
            neighbor_idx = neighbor_atom.GetIdx()
            if neighbor_atom.GetAtomicNum() == 6:
                # Check if neighbor atom is part of sugar ring
                for sugar in sugar_matches:
                    if neighbor_idx in sugar:
                        # Found potential glycosidic bond
                        # Now check for deoxyribose (no OH at 2' position)
                        # Identify the carbons in the sugar ring
                        carbon_indices = [idx for idx in sugar if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
                        # Check for OH group at 2' position
                        has_2prime_OH = False
                        for c_idx in carbon_indices:
                            c_atom = mol.GetAtomWithIdx(c_idx)
                            # Count how many carbons in ring the carbon atom is connected to
                            ring_carbons = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetIdx() in carbon_indices]
                            if len(ring_carbons) == 2:
                                # This is likely C2' position
                                for neighbor in c_atom.GetNeighbors():
                                    if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in sugar:
                                        # Found OH group at 2' position
                                        has_2prime_OH = True
                        if has_2prime_OH:
                            return False, "Sugar is ribose (has OH at 2' position)"
                        else:
                            return True, "Contains pyrimidine base attached to deoxyribose via N-glycosidic bond"

    return False, "No glycosidic bond between pyrimidine base and sugar found"