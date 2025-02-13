"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is defined as a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved biphenyl pattern with two connected benzene rings
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    matches = mol.GetSubstructMatches(biphenyl_pattern)

    if not matches:
        return False, "No biphenyl structure found"

    # Check each biphenyl match and count chlorines specifically on these phenyl rings
    chlorine_count = 0
    chlorine_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17}

    for match in matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'C' and atom.IsInRing():
                # Check neighbors for attached chlorines
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in chlorine_atoms:
                        chlorine_count += 1

    # Check if the chlorine count is between 2 and 10
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Contains {chlorine_count} chlorine atoms on the biphenyl structure, needs between 2 and 10"

    return True, "Contains biphenyl structure with acceptable number of chlorine atoms (between 2 and 10)"