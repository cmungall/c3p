"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,O])([O-,O])[O-,O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for serine moiety (NH2-CH-COOH)
    serine_pattern = Chem.MolFromSmarts("[NH2,NH3+][CH]C(=O)[OH,O-]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine moiety found"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Basic element count checks
    if c_count < 10:
        return False, "Too few carbons for phosphatidyl-L-serine"
    if o_count < 8:
        return False, "Too few oxygens for phosphatidyl-L-serine"
    if n_count != 1:
        return False, "Must have exactly one nitrogen (from serine)"
    if p_count != 1:
        return False, "Must have exactly one phosphorus"

    # Check for phosphoserine linkage
    phosphoserine_pattern = Chem.MolFromSmarts("[P](=O)([O-,O])([O-,O])OC[C@H](N)C(=O)[O-,OH]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Missing characteristic phosphoserine linkage"

    return True, "Contains phosphatidyl group esterified to serine with characteristic fatty acid chains"