"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:51753 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens are replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitro groups (-NO2)
    # More specific SMARTS pattern for nitro groups
    nitro_pattern = Chem.MolFromSmarts('[NX3+](=[OX1])([O-])-[#6]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Get all atoms
    atoms = mol.GetAtoms()
    
    # Check if molecule contains only C, H, N, and O
    allowed_atoms = {1, 6, 7, 8}  # H, C, N, O atomic numbers
    for atom in atoms:
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains elements other than C, H, N, O"

    # Check for carbon atoms (must have at least one)
    carbon_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Check for unwanted functional groups that would make it not a pure nitrohydrocarbon
    unwanted_groups = [
        ('[CX3](=O)[OX2H1]', 'carboxylic acid'),  # More specific SMARTS
        ('[CX3](=O)[O-]', 'carboxylate'),
        ('[CX3]=O', 'ketone/aldehyde'),
        ('[NX3;H2,H1;!$(NC=O)]', 'amine'),  # More specific, excludes nitro
        ('[SX2]', 'thiol/sulfide'),
        ('[Cl]', 'chloride'),
        ('[Br]', 'bromide'),
        ('[F]', 'fluoride'),
        ('[I]', 'iodide'),
        ('[OX2H]', 'hydroxyl'),  # More specific SMARTS
        ('[OX2H]-[CX4]', 'alcohol'),  # More specific SMARTS
    ]
    
    for pattern, group_name in unwanted_groups:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
            return False, f"Contains {group_name} group"

    # Count nitro groups
    num_nitro = len(nitro_matches)
    
    # Verify all nitrogen atoms are part of nitro groups
    n_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 7]
    if len(n_atoms) != len(nitro_matches):
        return False, "Contains nitrogen atoms not in nitro groups"
    
    # Verify all oxygen atoms are part of nitro groups
    o_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 8]
    if len(o_atoms) != 2 * len(nitro_matches):
        return False, "Contains oxygen atoms not in nitro groups"
    
    return True, f"Hydrocarbon with {num_nitro} nitro group{'s' if num_nitro > 1 else ''}"