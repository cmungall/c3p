"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:18035 diradylglycerol
A diradylglycerol is defined as any lipid that is glycerol bearing two substituent groups - 
either acyl, alkyl, or alk-1-enyl - at any two of the three possible positions.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.

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

    # Look for glycerol backbone pattern (C-C-C with 2 or 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X3,CH2X4][CHX4][CH2X3,CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 2 substituents (acyl, alkyl, or alk-1-enyl)
    substituent_pattern = Chem.MolFromSmarts("[OX2][CX3]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if len(substituent_matches) != 2:
        return False, f"Found {len(substituent_matches)} substituents, need exactly 2"

    # Check for alkyl/alkenyl chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, f"Missing alkyl/alkenyl chains, got {len(chain_matches)}"

    return True, "Contains glycerol backbone with 2 substituent groups (acyl, alkyl, or alk-1-enyl)"