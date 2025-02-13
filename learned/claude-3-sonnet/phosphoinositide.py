"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:18035 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the
    hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Look for inositol ring attached to glycerol via phosphate
    inositol_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring attached to glycerol via phosphate"
    
    # Check for additional phosphate groups on inositol
    inositol_p_pattern = Chem.MolFromSmarts("OC1C(O)C(OP)C(O)C(O)C1OP")
    if not mol.HasSubstructMatch(inositol_p_pattern):
        return False, "Inositol ring not phosphorylated"

    return True, "Contains glycerol backbone with 2 fatty acid chains, inositol ring attached via phosphate, and additional phosphate groups on inositol"