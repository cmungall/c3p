"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: Beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is a lactam in which the amide bond is contained within a four-membered ring,
    which includes the amide nitrogen and the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over all rings
    for ring in atom_rings:
        # Check for four-membered rings
        if len(ring) == 4:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Count nitrogen atoms in the ring
            n_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 7]
            # Count carbon atoms in the ring
            c_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
            # Check if there is exactly one nitrogen and three carbons
            if len(n_atoms) == 1 and len(c_atoms) == 3:
                nitrogen = n_atoms[0]
                # Look for a carbonyl carbon connected to the nitrogen within the ring
                for neighbor in nitrogen.GetNeighbors():
                    if neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 6:
                        # Check if this carbon has a double bond to oxygen
                        for bond in neighbor.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                other_atom = bond.GetOtherAtom(neighbor)
                                if other_atom.GetAtomicNum() == 8:
                                    return True, "Molecule contains a beta-lactam ring"
    return False, "No beta-lactam ring found"