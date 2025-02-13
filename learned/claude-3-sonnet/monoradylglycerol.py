"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
A lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic property checks
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight < 100 or mol_weight > 500:
        return False, f"Molecular weight {mol_weight:.1f} outside expected range (100-500)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:
        return False, f"Too few carbons ({c_count}) for monoradylglycerol"
    if not (3 <= o_count <= 4):
        return False, f"Invalid number of oxygen atoms ({o_count}), expected 3-4"

    # Look for glycerol backbone - more general pattern
    glycerol_pattern = Chem.MolFromSmarts("[OX2,OH1][CH2][CH]([OX2,OH1])[CH2][OX2,OH1]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count free hydroxyls
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches != 2:
        return False, f"Found {oh_matches} free hydroxyl groups, need exactly 2"

    # Define substituent patterns
    patterns = {
        'acyl': [
            # 1/3-acyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CX3](=[OX1])[#6]",
            # 2-acyl
            "[OH1][CH2][CH]([OX2][CX3](=[OX1])[#6])[CH2][OH1]"
        ],
        'alkyl': [
            # 1/3-alkyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CX4;!$(C=O)]",
            # 2-alkyl
            "[OH1][CH2][CH]([OX2][CX4;!$(C=O)])[CH2][OH1]"
        ],
        'alkenyl': [
            # 1/3-alkenyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CH2][CH]=[CH]",
            # 2-alkenyl
            "[OH1][CH2][CH]([OX2][CH2][CH]=[CH])[CH2][OH1]"
        ]
    }

    # Count substituents
    total_matches = 0
    substituent_type = None
    
    for subst_type, subst_patterns in patterns.items():
        for pattern in subst_patterns:
            matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
            if matches > 0:
                total_matches += matches
                substituent_type = subst_type

    if total_matches != 1:
        return False, f"Found {total_matches} substituents, need exactly 1"

    # Additional validation for chain length
    chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    chain_matches = len(mol.GetSubstructMatches(chain_pattern))
    if chain_matches < 1:
        return False, "Substituent chain too short"

    return True, f"Valid monoradylglycerol with one {substituent_type} substituent"