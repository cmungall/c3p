"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: phosphatidylinositol compounds
Definition: Any glycerophosphoinositol having one phosphatidyl group esterified to one of the hydroxy groups of inositol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has:
    - A myo-inositol ring (cyclohexane with 6 OH groups)
    - Connected via phosphate to a glycerol backbone
    - Two fatty acid chains attached to the glycerol
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple SMARTS patterns for inositol ring to handle different representations
    inositol_patterns = [
        # Basic cyclohexane with 6 oxygens
        "[C]1[C][C][C][C][C]1([O,OH])[O,OH]([O,OH])[O,OH]([O,OH])[O,OH]",
        # More specific myo-inositol pattern
        "[CH1,CH2]1[CH1,CH2]([OH1,OH0])[CH1,CH2]([OH1,OH0])[CH1,CH2]([OH1,OH0])[CH1,CH2]([OH1,OH0])[CH1,CH2]([OH1,OH0])O1",
        # Pattern with explicit stereochemistry
        "[C@@H,C@H,CH]1[C@@H,C@H,CH]([OH])[C@@H,C@H,CH]([OH])[C@@H,C@H,CH]([OH])[C@@H,C@H,CH]([OH])[C@@H,C@H,CH]1[OH]"
    ]
    
    inositol_found = False
    for pattern in inositol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            inositol_found = True
            break
    
    if not inositol_found:
        return False, "No inositol ring found"

    # Look for phosphate group connected to inositol and glycerol
    phosphate_patterns = [
        "O[P](=O)(O)OC",
        "[OH0,OH1][P](=[O])(O[C])[O][C]",
        "[O,OH][P]([O,OH])(=O)OC"
    ]
    
    phosphate_found = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            phosphate_found = True
            break
            
    if not phosphate_found:
        return False, "No phosphate group found"

    # Look for glycerol backbone with flexible pattern
    glycerol_patterns = [
        "[CH2X4,CH2][CHX4,CH][CH2X4,CH2]",
        "CCC([OH1,OH0,O])",
        "[CH2][CH]([O])[CH2]"
    ]
    
    glycerol_found = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            glycerol_found = True
            break
            
    if not glycerol_found:
        return False, "No glycerol backbone found"

    # Look for two ester groups
    ester_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains
    fatty_chain_patterns = [
        "[CH2][CH2][CH2][CH2]",  # saturated
        "[CH2][CH]=[CH][CH2]",   # unsaturated
        "CCCC"                    # generic
    ]
    
    chain_count = 0
    for pattern in fatty_chain_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        chain_count += len(matches)
        if chain_count >= 2:
            break
    
    if chain_count < 2:
        return False, "Missing fatty acid chains"

    # Verify overall composition
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 1:
        return False, f"Should have exactly 1 phosphorus atom, found {phosphorus_count}"

    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 11:
        return False, f"Too few oxygen atoms for phosphatidylinositol (found {oxygen_count})"

    return True, "Contains inositol ring connected via phosphate to glycerol backbone with two fatty acid chains"