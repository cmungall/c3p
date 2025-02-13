"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: monoacylglycerol
A glyceride in which any one of the R groups (position not specified) is an acyl group 
while the remaining two R groups can be either H or alkyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with one fatty acid chain attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate groups - exclude phospholipids
    phosphate_pattern = Chem.MolFromSmarts("[PX4]")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group - not a monoacylglycerol"

    # More comprehensive glycerol backbone patterns
    glycerol_patterns = [
        # Basic pattern with any position ester
        "[OX2,OH1]-[CH2X4]-[CHX4]-[CH2X4]-[OX2,OH1]",
        # 1- or 3-position ester
        "[CX3](=[OX1])[OX2]-[CH2X4]-[CHX4](-[OX2,OH1])-[CH2X4]-[OX2,OH1]",
        # 2-position ester
        "[OX2,OH1]-[CH2X4]-[CHX4](-[OX2][CX3]=[OX1])-[CH2X4]-[OX2,OH1]"
    ]
    
    has_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_glycerol = True
            break
            
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count free or alkylated hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H1,OX2][CH2X4,CHX4]-[CHX4]-[CH2X4]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) < 2:
        return False, "Insufficient free or substituted hydroxyl groups"

    # Verify no other complex modifications
    complex_pattern = Chem.MolFromSmarts("[N,S,P,F,Cl,Br,I]")
    if mol.HasSubstructMatch(complex_pattern):
        return False, "Contains non-standard atoms or groups"

    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:
        return False, "Too few carbons for monoacylglycerol"
    if o_count != 4:  # Must have exactly 4 oxygens
        return False, "Must have exactly 4 oxygens (1 ester + 2 hydroxyls)"

    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120 or mol_wt > 1000:
        return False, "Molecular weight outside reasonable range for monoacylglycerol"

    # Verify acyl chain connection - multiple patterns to catch all positions
    acyl_patterns = [
        # 1- or 3-position
        "[CX3](=[OX1])[OX2]-[CH2X4]-[CHX4](-[OX2,OH1])-[CH2X4]-[OX2,OH1]",
        # 2-position
        "[OX2,OH1]-[CH2X4]-[CHX4](-[OX2][CX3]=[OX1])-[CH2X4]-[OX2,OH1]",
        # Alternative patterns
        "[CX3](=[OX1])[OX2]-[CH2X4,CHX4]-[CH1X4](-[OX2,OH1])-[CH2X4]-[OX2,OH1]",
        "[OX2,OH1]-[CH2X4]-[CH1X4](-[OX2][CX3]=[OX1])-[CH2X4]-[OX2,OH1]"
    ]
    
    for pattern in acyl_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, "Contains glycerol backbone with one acyl group attached via ester bond"

    return False, "Acyl chain not properly connected to glycerol backbone"