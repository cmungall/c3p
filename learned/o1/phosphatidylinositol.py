"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: Phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_inositol_ring(mol):
    """
    Checks if the molecule contains an inositol ring (a cyclohexane ring with six hydroxyl groups).

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        set: Set of atom indices forming the inositol ring, or empty set if not found
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for ring in atom_rings:
        if len(ring) == 6:
            is_inositol = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    is_inositol = False
                    break
                has_oxygen = False
                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                        has_oxygen = True
                        break
                if not has_oxygen:
                    is_inositol = False
                    break
            if is_inositol:
                return set(ring)  # Return set of atom indices of the inositol ring
    return set()

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol having one phosphatidyl group
    esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import ChemicalFeatures
    from rdkit.Chem.Pharm2D import Generate

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for two ester groups (C(=O)O-C)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, expected at least 2"

    # Check for inositol ring
    inositol_atoms = has_inositol_ring(mol)
    if not inositol_atoms:
        return False, "No inositol ring found"

    # Check for phosphate group (P(=O)(O)(O))
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check if phosphate is connected to inositol ring
    phosphate_connected_to_inositol = False
    for match in phosphate_matches:
        p_idx = match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() == 8:
                for neighbor in o_atom.GetNeighbors():
                    if neighbor.GetIdx() in inositol_atoms:
                        phosphate_connected_to_inositol = True
                        break
            if phosphate_connected_to_inositol:
                break
        if phosphate_connected_to_inositol:
            break
    if not phosphate_connected_to_inositol:
        return False, "Phosphate group not connected to inositol ring"

    # Check for glycerol backbone connected to phosphate group
    glycerol_pattern = Chem.MolFromSmarts("COC(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check if phosphate is connected to glycerol backbone
    phosphate_connected_to_glycerol = False
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    for match in phosphate_matches:
        p_idx = match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() == 8:
                for neighbor in o_atom.GetNeighbors():
                    if neighbor.GetIdx() in glycerol_atoms:
                        phosphate_connected_to_glycerol = True
                        break
            if phosphate_connected_to_glycerol:
                break
        if phosphate_connected_to_glycerol:
            break
    if not phosphate_connected_to_glycerol:
        return False, "Phosphate group not connected to glycerol backbone"

    # All checks passed
    return True, "Molecule is a phosphatidylinositol"

# Example usage:
# smiles = "[C@@H]1(C(C(C([C@H](C1O)O)O)O)O)OP(OC[C@@](COC(CCCCCCCCCCCCCCC)=O)(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)[H])(O)=O"
# result, reason = is_phosphatidylinositol(smiles)
# print(result, reason)