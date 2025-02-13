"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    PIPs contain a myo-inositol ring with at least one phosphate group attached directly
    to the ring hydroxyl groups (not counting the phosphate that links to glycerol).
    Most common forms are PI(3)P, PI(4)P, PI(5)P, PI(3,4)P2, PI(3,5)P2, PI(4,5)P2, and PI(3,4,5)P3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring (6-membered ring with all carbons and correct substitution)
    inositol_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Find inositol ring atoms
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"
    inositol_atoms = set(inositol_matches[0])

    # Pattern for phosphate group attached to oxygen
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"

    # Find all phosphates attached to inositol ring
    inositol_phosphates = 0
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    for match in phosphate_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        # Get the atom the oxygen is connected to
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIdx() in inositol_atoms:
                inositol_phosphates += 1
                break

    # Must have at least one phosphate directly on inositol ring
    # (not counting the bridging phosphate to glycerol)
    if inositol_phosphates < 2:  # One phosphate is the bridge to glycerol
        return False, "No additional phosphate groups on inositol ring beyond the glycerol linkage"

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    # Verify either glycerol backbone or modified structure
    is_modified = False
    if not mol.HasSubstructMatch(glycerol_pattern):
        is_modified = True
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "Neither glycerol backbone nor modified structure found"

    # For standard PIPs, verify fatty acid chains
    if not is_modified:
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
        if fatty_acid_matches < 2:
            return False, "Missing fatty acid chains in glycerol-based PIP"

    # Classify based on number of phosphates on inositol
    actual_modifications = inositol_phosphates - 1  # Subtract the bridging phosphate
    
    if actual_modifications >= 3:
        return True, "Contains inositol ring with trisphosphate modification (PIP3)"
    elif actual_modifications == 2:
        return True, "Contains inositol ring with bisphosphate modification (PIP2)"
    elif actual_modifications == 1:
        return True, "Contains inositol ring with monophosphate modification (PIP)"
    else:
        return False, "Insufficient phosphate modifications on inositol ring"