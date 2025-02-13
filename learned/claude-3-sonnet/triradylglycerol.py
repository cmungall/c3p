"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:35781 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol backbone with three substituents - either acyl, alkyl, or alk-1-enyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for 3 substituents (acyl, alkyl, or alk-1-enyl groups)
    acyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    alkyl_pattern = Chem.MolFromSmarts("[CX4]")
    alkenyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    
    substituents = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetTotalNumHs() == 0:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.HasSubstructMatch(acyl_pattern):
                substituents.append("acyl")
            elif neighbor.HasSubstructMatch(alkyl_pattern):
                substituents.append("alkyl")
            elif neighbor.HasSubstructMatch(alkenyl_pattern):
                substituents.append("alk-1-enyl")
    
    if len(set(substituents)) != 3 or len(substituents) != 3:
        return False, f"Found {len(substituents)} substituents, need exactly 3 (acyl, alkyl, alk-1-enyl)"
    
    return True, "Glycerol backbone with 3 substituents: acyl, alkyl, and alk-1-enyl"