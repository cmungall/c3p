"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
"""
from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is any glycerophospholipid having the polar alcohol inositol
    esterified to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find phosphorus atoms in the molecule
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphate group found"

    # Function to check if an atom is part of an inositol ring
    def is_inositol_ring(atom):
        # Inositol ring: 6-membered carbocycle with hydroxyl groups
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) != 6:
                continue
            if atom.GetIdx() not in ring:
                continue
            # Check if all atoms in ring are carbon
            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                continue
            # Check for hydroxyl groups on ring carbons
            hydroxyl_count = 0
            for idx in ring:
                ring_atom = mol.GetAtomWithIdx(idx)
                has_oh = False
                for neighbor in ring_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        has_oh = True
                        break
                if has_oh:
                    hydroxyl_count += 1
            if hydroxyl_count >= 5:
                return True
        return False

    # Function to check if an atom is part of a glycerol backbone
    def is_glycerol_backbone(atom):
        # Glycerol backbone: 3 connected carbons with specific attachments
        glycerol_smarts = "[CX4H2][CX4H](O)[CX4H2](O)"
        glycerol = Chem.MolFromSmarts(glycerol_smarts)
        matches = mol.GetSubstructMatches(glycerol)
        return any(atom.GetIdx() in match for match in matches)

    # Check each phosphorus atom
    for p_atom in p_atoms:
        o_neighbors = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        connected_to_inositol = False
        connected_to_glycerol = False

        for o_atom in o_neighbors:
            # Skip P=O double bond
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_atom.GetIdx())
            if bond.GetBondTypeAsDouble() == 2.0:
                continue  # This is a double bond (P=O)
            # Check the other atom connected to this oxygen
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetIdx() == p_atom.GetIdx():
                    continue
                if is_inositol_ring(neighbor):
                    connected_to_inositol = True
                elif is_glycerol_backbone(neighbor):
                    connected_to_glycerol = True

        if connected_to_inositol and connected_to_glycerol:
            return True, "Molecule is a glycerophosphoinositol"

    return False, "Phosphate group not connected to both inositol and glycerol backbone"

__metadata__ = {
    'chemical_class': {
        'name': 'glycerophosphoinositol',
        'definition': 'Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.'
    }
}