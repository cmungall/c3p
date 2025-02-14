"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:24048 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups - either acyl, alkyl, or alk-1-enyl - at any two of the three possible positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]([OX2])[CHX4]([OX2])[CH2X4]([OX2])")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No glycerol backbone found"

    glycerol_atoms = glycerol_match

    # Count attachments
    attachment_points = 0

    # Ester attachment pattern, explicitly attached to glycerol carbon
    ester_attach_pattern = Chem.MolFromSmarts("[CH2X4;!H0,CHX4;!H0]-[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_attach_pattern)
    attachment_points += len(ester_matches)


    # Ether attachment (alkyl) pattern, explicitly attached to glycerol carbon
    ether_attach_pattern = Chem.MolFromSmarts("[CH2X4;!H0,CHX4;!H0]-[OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_attach_pattern)
    attachment_points += len(ether_matches)

    # Enol ether attachment (alk-1-enyl), explicitly attached to glycerol carbon
    enol_ether_attach_pattern = Chem.MolFromSmarts("[CH2X4;!H0,CHX4;!H0]-[OX2][CX3]=[CX3]")
    enol_ether_matches = mol.GetSubstructMatches(enol_ether_attach_pattern)
    attachment_points += len(enol_ether_matches)

    if attachment_points != 2:
         return False, f"Found {attachment_points} substituents, need exactly 2"
    
    return True, "Contains glycerol backbone with exactly two substituents (acyl, alkyl or alk-1-enyl) attached"