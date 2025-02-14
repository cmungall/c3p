"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: aldose
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is an aldehydic parent sugar (polyhydroxy aldehyde) or its intramolecular hemiacetal form
    (a cyclic structure like a furanose or pyranose ring).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule contains only C, H, and O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1, 8):
            return False, "Contains elements other than C, H, and O"

    # Check number of carbon atoms (aldoses typically have 3 or more carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Too few carbons ({c_count}) for an aldose"

    # Define patterns
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")  # Aldehyde group
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")  # Ketone group (exclude ketoses)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")  # Carboxylic acid group
    hemiacetal_oxygen_pattern = Chem.MolFromSmarts("[OX2H][CX4H]")  # Hydroxyl oxygen attached to carbon

    # Exclude molecules with ketone groups (ketoses) or carboxylic acids (uronic acids)
    if mol.HasSubstructMatch(ketone_pattern):
        return False, "Contains ketone group, likely a ketose"
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains carboxylic acid group, likely a uronic acid or aldonic acid"

    # Check for aldehyde group (open-chain form)
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Open-chain aldose
        # Count hydroxyl groups attached to carbons (excluding aldehyde carbon)
        hydroxyl_pattern = Chem.MolFromSmarts("[CX4;$([CH0,CH1,CH2][OH])]")  # Carbon with hydroxyl
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        num_hydroxyls = len(hydroxyl_matches)

        # Require at least two hydroxyl groups for aldoses with three carbons
        min_hydroxyls = max(2, c_count - 2)  # Allow for deoxy sugars
        if num_hydroxyls < min_hydroxyls:
            return False, f"Too few hydroxyl groups ({num_hydroxyls}) for an aldose"

        return True, "Open-chain aldose with aldehyde group and sufficient hydroxylated carbons"

    else:
        # Check for cyclic hemiacetal forms (furanose or pyranose rings)
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 0:
            return False, "Does not contain aldehyde group or cyclic hemiacetal ring"

        # Find rings of size 5 or 6 containing one oxygen (hemiacetal form)
        found_ring = False
        rings = ring_info.AtomRings()
        for ring in rings:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_size = len(ring)
            if ring_size not in (5, 6):
                continue
            o_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
            if len(o_atoms) != 1:
                continue
            o_atom = o_atoms[0]
            # Check if oxygen is connected to anomeric carbon
            o_neighbors = o_atom.GetNeighbors()
            if len(o_neighbors) != 2:
                continue
            # Anomeric carbon is connected to two oxygens
            for nbr in o_neighbors:
                if nbr.GetAtomicNum() == 6:
                    o_nb = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 8 and a.GetIdx() != o_atom.GetIdx()]
                    if len(o_nb) >= 1:
                        anomeric_c = nbr
                        found_ring = True
                        break
            if found_ring:
                break

        if not found_ring:
            return False, "No suitable furanose or pyranose ring with hemiacetal found"

        # Ensure no ketone groups are present
        if mol.HasSubstructMatch(ketone_pattern):
            return False, "Contains ketone group, likely a ketose"

        # Exclude molecules with carboxylic acids
        if mol.HasSubstructMatch(carboxylic_acid_pattern):
            return False, "Contains carboxylic acid group, likely a uronic acid"

        # Count hydroxyl groups attached to carbons (excluding ring oxygen)
        hydroxyl_pattern = Chem.MolFromSmarts("[CX4H][OH]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        num_hydroxyls = len(hydroxyl_matches)

        # Require at least two hydroxyl groups for cyclic aldoses
        min_hydroxyls = max(2, c_count - 3)  # Allow for deoxy sugars and ring closure
        if num_hydroxyls < min_hydroxyls:
            return False, f"Too few hydroxyl groups ({num_hydroxyls}) for a cyclic aldose"

        return True, f"Cyclic aldose (hemiacetal form) with ring size {ring_size} and sufficient hydroxyl groups"

    return False, "Does not match aldose criteria"