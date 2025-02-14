"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound,
    meaning the nitrile group (-C#N) is attached to an aliphatic (sp3-hybridized, non-aromatic) carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile group SMARTS pattern (carbon triple-bonded to nitrogen)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    # Check for presence of nitrile groups
    if not nitrile_matches:
        return False, "No nitrile group found"

    # Flag to track if any nitrile is attached to an aliphatic carbon
    has_aliphatic_nitrile = False

    # Iterate over all nitrile groups in the molecule
    for match in nitrile_matches:
        nitrile_c_idx = match[0]  # Index of carbon in nitrile group
        nitrile_n_idx = match[1]  # Index of nitrogen in nitrile group
        nitrile_c = mol.GetAtomWithIdx(nitrile_c_idx)
        nitrile_n = mol.GetAtomWithIdx(nitrile_n_idx)

        # Get neighbors of the nitrile carbon atom (exclude the nitrogen)
        neighbors = [atom for atom in nitrile_c.GetNeighbors() if atom.GetIdx() != nitrile_n_idx]

        # If nitrile carbon has no other neighbors, it's a terminal atom (unlikely but check)
        if not neighbors:
            continue  # No adjacent atom to nitrile carbon other than nitrogen

        # Check if the adjacent atom to nitrile carbon is:
        # - A carbon atom
        # - Not aromatic
        # - sp3 hybridized
        # - Not double or triple bonded to nitrile carbon
        for neighbor in neighbors:
            # Check if neighbor is a carbon atom
            if neighbor.GetAtomicNum() != 6:
                continue  # Not a carbon atom

            # Check if neighbor is aromatic
            if neighbor.GetIsAromatic():
                continue  # Neighbor is aromatic

            # Check hybridization of neighbor (should be sp3)
            if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue  # Neighbor is not sp3 hybridized

            # Check bond type between nitrile carbon and neighbor
            bond = mol.GetBondBetweenAtoms(nitrile_c_idx, neighbor.GetIdx())
            if bond is None:
                continue  # No bond found (should not happen)

            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue  # Bond is not a single bond

            # Passed all checks; nitrile is attached to an aliphatic carbon
            has_aliphatic_nitrile = True
            break  # No need to check other neighbors

        if has_aliphatic_nitrile:
            break  # Found an aliphatic nitrile group

    if has_aliphatic_nitrile:
        return True, "Contains nitrile group(s) attached to aliphatic carbon(s)"
    else:
        return False, "No nitrile groups attached to aliphatic carbons found"