"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is any aldoxime derived from an aliphatic aldehyde,
    meaning it contains an aldoxime group (R-CH=N-OH) attached to an aliphatic chain,
    with no aromatic rings or cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic atoms
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms"

    # Check for any rings (acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures"

    # Define more specific aldoxime SMARTS pattern (R-CH=N-OH)
    aldoxime_pattern = Chem.MolFromSmarts('[CH1]=[N][OH1]')
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime functional group found"

    # Iterate over matches to ensure no extra substituents on N or O
    for match in matches:
        carbon = mol.GetAtomWithIdx(match[0])
        nitrogen = mol.GetAtomWithIdx(match[1])
        oxygen = mol.GetAtomWithIdx(match[2])

        # Check that nitrogen has only two neighbors (carbon and oxygen)
        if nitrogen.GetDegree() != 2:
            continue

        # Check that oxygen has only one neighbor (nitrogen)
        if oxygen.GetDegree() != 1:
            continue

        # If all checks pass, return True
        return True, "Contains aldoxime group derived from an aliphatic aldehyde"

    return False, "No valid aldoxime group found"