"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:35741 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol compound having one of three possible substituent groups
    (acyl, alkyl, or alk-1-enyl) at each of the three possible positions sn-1, sn-2, or sn-3.

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

    # Look for 3 ester groups (-O-C(=O)-), alkyl chains (-C-), or alkenyl chains (-C=C-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    alkyl_pattern = Chem.MolFromSmarts("[CX4]~[CX4]")
    alkenyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    
    substituents = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O":
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, atom.GetIdx())
            if env.HasSubstructMatch(ester_pattern):
                substituents.append("acyl")
            elif env.HasSubstructMatch(alkyl_pattern):
                substituents.append("alkyl")
            elif env.HasSubstructMatch(alkenyl_pattern):
                substituents.append("alk-1-enyl")

    if len(substituents) != 3:
        return False, f"Found {len(substituents)} substituents, need exactly 3"
    
    # Check for valid substituent combinations
    valid_combos = [
        ["acyl", "acyl", "acyl"],
        ["acyl", "acyl", "alkyl"],
        ["acyl", "acyl", "alk-1-enyl"],
        ["acyl", "alkyl", "alkyl"],
        ["acyl", "alkyl", "alk-1-enyl"],
        ["acyl", "alk-1-enyl", "alk-1-enyl"],
        ["alkyl", "alkyl", "alkyl"],
        ["alkyl", "alkyl", "alk-1-enyl"],
        ["alkyl", "alk-1-enyl", "alk-1-enyl"],
        ["alk-1-enyl", "alk-1-enyl", "alk-1-enyl"]
    ]
    
    if sorted(substituents) not in valid_combos:
        return False, f"Invalid substituent combination: {', '.join(substituents)}"

    return True, "Contains glycerol backbone with 3 substituents: acyl, alkyl, or alk-1-enyl"