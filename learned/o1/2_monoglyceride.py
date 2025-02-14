"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone where the middle carbon (position 2)
    is esterified with a fatty acid, and the other two carbons (positions 1 and 3)
    have hydroxyl groups and no other substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern: three connected carbons
    glycerol_pattern = Chem.MolFromSmarts("C-C-C")
    matches = mol.GetSubstructMatches(glycerol_pattern)

    if not matches:
        return False, "No glycerol backbone found"

    # Iterate over glycerol backbone matches
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check that c1 and c3 have hydroxyl groups and no other substituents
        def is_terminal_carbon_ok(carbon_atom, connected_idx):
            oh_found = False
            for neighbor in carbon_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == connected_idx:
                    continue  # Skip the bond to the central carbon
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    # Oxygen with only one bond (hydroxyl group)
                    oh_found = True
                else:
                    return False  # Other substituent found
            return oh_found

        c1_ok = is_terminal_carbon_ok(c1, c2_idx)
        c3_ok = is_terminal_carbon_ok(c3, c2_idx)

        if not (c1_ok and c3_ok):
            continue  # Positions 1 and 3 must have only hydroxyl groups

        # Check that c2 is connected to an ester group and no other substituents
        ester_found = False
        for neighbor in c2.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in [c1_idx, c3_idx]:
                continue  # Skip bonds to the glycerol backbone carbons
            if neighbor.GetAtomicNum() == 8:
                # Check if this oxygen is part of an ester linkage
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetIdx() == c2_idx:
                        continue  # Skip the bond back to c2
                    if neighbor2.GetAtomicNum() == 6:
                        # Check for carbonyl (C=O)
                        for bond in neighbor2.GetBonds():
                            other_atom = bond.GetOtherAtom(neighbor2)
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and other_atom.GetAtomicNum() == 8:
                                ester_found = True
                                break
                        if ester_found:
                            break
                if ester_found:
                    continue
                else:
                    return False, "Position 2 oxygen is not part of ester linkage"
            else:
                return False, "Position 2 has unexpected substituents"

        if not ester_found:
            continue  # Ester linkage at position 2 is required

        # Check that there is only one ester group in the molecule
        ester_pattern = Chem.MolFromSmarts("C(=O)O")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) != 1:
            return False, f"Found {len(ester_matches)} ester groups, expected 1"

        # Ensure there are exactly two hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) != 2:
            return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, expected 2"

        # Passed all checks
        return True, "Molecule is a 2-monoglyceride with acyl group at position 2"

    return False, "Molecule does not match the 2-monoglyceride structure"