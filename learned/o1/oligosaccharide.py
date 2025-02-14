"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound in which monosaccharide units are joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for monosaccharide rings (5 or 6 membered ring with oxygen)
    monosaccharide_smarts = "[C,O]1[C,O][C,O][C,O][C,O][C,O]1"  # Simplified pattern
    monosaccharide_pattern = Chem.MolFromSmarts(monosaccharide_smarts)

    # Find all matches for monosaccharide units
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    if not monosaccharide_matches:
        return False, "No monosaccharide units found"

    monosaccharide_atoms = [set(match) for match in monosaccharide_matches]
    num_monosaccharides = len(monosaccharide_atoms)

    # Build atom index to monosaccharide index mapping
    atom_mono_map = {}
    for idx, atom_set in enumerate(monosaccharide_atoms):
        for atom_idx in atom_set:
            atom_mono_map.setdefault(atom_idx, []).append(idx)

    # Detect glycosidic linkages between monosaccharide units
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()

        # Skip if both atoms are in the same monosaccharide unit
        mono_indices1 = set(atom_mono_map.get(idx1, []))
        mono_indices2 = set(atom_mono_map.get(idx2, []))
        if not mono_indices1 or not mono_indices2:
            continue
        if mono_indices1 == mono_indices2:
            continue

        # Check for oxygen bridges between monosaccharides (glycosidic linkages)
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                glycosidic_bonds.append(bond)

    num_glycosidic_bonds = len(glycosidic_bonds)
    if num_glycosidic_bonds < num_monosaccharides - 1:
        return False, "Not all monosaccharide units are connected via glycosidic linkages"

    # Ensure the number of monosaccharide units is within the typical range for oligosaccharides
    if num_monosaccharides < 2:
        return False, "Less than 2 monosaccharide units found"
    elif num_monosaccharides > 20:
        return False, "More than 20 monosaccharide units found; may be a polysaccharide"

    return True, f"Contains {num_monosaccharides} monosaccharide units connected via glycosidic linkages"