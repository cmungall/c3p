"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to exactly two carbon atoms and one hydrogen,
    excluding cases where nitrogen is part of amides, aromatic systems, nitro groups,
    nitriles, imines, or other functional groups where nitrogen is double-bonded
    or bonded to heteroatoms other than hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Skip if not nitrogen
            continue

        # Check if nitrogen has exactly one hydrogen
        if atom.GetTotalNumHs(includeNeighbors=True) != 1:
            continue

        # Check if nitrogen is bonded to exactly two carbon atoms
        neighbors = atom.GetNeighbors()
        carbon_count = 0
        has_heteroatom_neighbor = False
        for nbr in neighbors:
            atomic_num = nbr.GetAtomicNum()
            if atomic_num == 6:
                carbon_count += 1
            elif atomic_num != 1:  # Exclude hydrogen
                has_heteroatom_neighbor = True
                break

        if carbon_count != 2 or has_heteroatom_neighbor:
            continue

        # Exclude nitrogens with double bonds (e.g., imines, nitriles)
        has_double_bond = False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetAtomicNum() != 1:
                has_double_bond = True
                break

        if has_double_bond:
            continue

        # Exclude nitrogens in amides (N-C(=O))
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Neighbor is carbon
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx()).GetBondTypeAsDouble() == 2.0:
                        is_amide = True
                        break
                if is_amide:
                    break
        if is_amide:
            continue

        # Exclude nitro groups (N(=O)-O)
        is_nitro = False
        if atom.GetFormalCharge() == 1:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    is_nitro = True
                    break
        if is_nitro:
            continue

        # Exclude aromatic nitrogens participating in aromatic systems
        if atom.GetIsAromatic():
            continue

        # If all checks passed, it's a secondary amine
        return True, "Contains secondary amine group"

    return False, "Does not contain secondary amine group"