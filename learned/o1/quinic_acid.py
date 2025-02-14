"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid - a cyclohexane ring with multiple hydroxyl groups (or their derivatives) and a carboxylic acid group attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cyclohexane ring
    cyclohexane_smarts = "[R]{6}"
    cyclohexane = Chem.MolFromSmarts(cyclohexane_smarts)
    ring_matches = mol.GetSubstructMatches(cyclohexane)
    if not ring_matches:
        return False, "No cyclohexane ring found"

    # For each cyclohexane ring found
    for match in ring_matches:
        ring_atoms = set(match)
        substituent_count = 0
        carboxylic_acid_found = False

        # Check substituents on ring atoms
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in ring_atoms:
                    # Check for carboxylic acid or ester group attached to ring
                    if neighbor.GetAtomicNum() == 6:
                        neighbor_atom = neighbor
                        # Look for -C(=O)O pattern
                        c_atom = neighbor_atom
                        o_double = False
                        o_single = False
                        for n in c_atom.GetNeighbors():
                            if n.GetIdx() != idx:
                                bond_type = c_atom.GetBondBetweenAtom(n.GetIdx()).GetBondType()
                                if n.GetAtomicNum() == 8:
                                    if bond_type == Chem.rdchem.BondType.DOUBLE:
                                        o_double = True
                                    elif bond_type == Chem.rdchem.BondType.SINGLE:
                                        o_single = True
                        if o_double and o_single:
                            carboxylic_acid_found = True
                            continue

                    # Check for hydroxyl group (–OH)
                    if neighbor.GetAtomicNum() == 8:
                        # Ensure it's a hydroxyl group and not an ether
                        if neighbor.GetDegree() == 1:
                            substituent_count += 1
                            continue
                        else:
                            # Check if oxygen is part of an ester linkage
                            o_atom = neighbor
                            for n in o_atom.GetNeighbors():
                                if n.GetIdx() != idx:
                                    if n.GetAtomicNum() == 6:
                                        # Possible ester linkage
                                        c_atom = n
                                        double_bonded_o = False
                                        for nn in c_atom.GetNeighbors():
                                            if nn.GetIdx() != o_atom.GetIdx():
                                                bond_type = c_atom.GetBondBetweenAtom(nn.GetIdx()).GetBondType()
                                                if nn.GetAtomicNum() == 8 and bond_type == Chem.rdchem.BondType.DOUBLE:
                                                    double_bonded_o = True
                                        if double_bonded_o:
                                            substituent_count +=1
                                            break

        if substituent_count >= 3 and carboxylic_acid_found:
            return True, "Contains cyclitol carboxylic acid core (quinic acid or derivative)"
        else:
            if substituent_count < 3:
                return False, f"Found {substituent_count} hydroxyl/ester groups on ring, need at least 3"
            if not carboxylic_acid_found:
                return False, "No carboxylic acid group attached to ring"

    return False, "Does not contain the required cyclitol carboxylic acid core"