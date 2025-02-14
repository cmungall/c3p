"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a lactone (cyclic ester) having a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate if gamma-lactone is found
    found_gamma_lactone = False

    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if the ring is five-membered
        if len(ring) == 5:
            # Get atoms in the ring
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Iterate over atoms in the ring to find ester linkage
            for atom in ring_atoms:
                # Check if atom is a carbon atom
                if atom.GetAtomicNum() == 6:
                    has_carbonyl_oxygen = False  # Flag for carbonyl oxygen
                    has_ester_oxygen = False     # Flag for ester oxygen
                    # Examine bonds of the carbon atom
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        # Check if neighbor atom is in the same ring
                        if neighbor.GetIdx() in ring:
                            # Check if neighbor is an oxygen atom
                            if neighbor.GetAtomicNum() == 8:
                                # Check bond type (double or single)
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    has_carbonyl_oxygen = True
                                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                    has_ester_oxygen = True
                    # If both carbonyl and ester oxygens are present, gamma-lactone is found
                    if has_carbonyl_oxygen and has_ester_oxygen:
                        found_gamma_lactone = True
                        return True, "Contains gamma-lactone ring (five-membered lactone ring)"

    # If no gamma-lactone ring is found after examining all rings
    return False, "No gamma-lactone ring found"