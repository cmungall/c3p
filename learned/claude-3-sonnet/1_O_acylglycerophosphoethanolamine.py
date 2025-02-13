"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    These compounds have:
    - A glycerol backbone
    - An acyl group (ester) at sn-1 position
    - A hydroxyl at sn-2 position
    - A phosphoethanolamine group at sn-3 position
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone connected to correct groups
    # [CH2X4]-[CHX4]-[CH2X4] where:
    # - One CH2 has ester (sn-1)
    # - CH has OH (sn-2)
    # - Other CH2 has phosphate (sn-3)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphoethanolamine group with primary amine
    # -P(=O)(O)-O-CH2-CH2-NH2/NH3+
    phosphoethanolamine = Chem.MolFromSmarts("[PX4](=O)([OH,O-])[OX2][CH2X4][CH2X4][NH2X3,NH3X4+]")
    if not mol.HasSubstructMatch(phosphoethanolamine):
        return False, "No phosphoethanolamine group with primary amine found"

    # Exclude phosphocholine (N(CH3)3+)
    phosphocholine = Chem.MolFromSmarts("[NX4+](C)(C)(C)")
    if mol.HasSubstructMatch(phosphocholine):
        return False, "Contains quaternary amine (phosphocholine) instead of primary amine"

    # Check for single ester group at sn-1
    # R-C(=O)-O-CH2- connected to glycerol
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH2X4][CHX4]([OX2H1])[CH2X4]O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group at sn-1 position found"
    
    # Count ester groups - should only be one
    ester_matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")))
    if ester_matches != 1:
        return False, f"Found {ester_matches} ester groups, should be exactly 1"

    # Exclude vinyl ethers (plasmalogens)
    vinyl_ether = Chem.MolFromSmarts("[OX2][CH2X4][CHX4]([OX2H1])[CH2X4]O[CH2X3]=[CH1X3]")
    if mol.HasSubstructMatch(vinyl_ether):
        return False, "Contains vinyl ether (plasmalogen) instead of ester at sn-1"

    # Exclude N-acylated versions
    n_acyl = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
    if mol.HasSubstructMatch(n_acyl):
        return False, "Contains N-acylation"

    # Verify basic composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if p_count != 1:
        return False, f"Should have exactly 1 phosphorus, found {p_count}"
    if o_count < 5:
        return False, f"Should have at least 5 oxygens, found {o_count}"
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Check for acyl chain length - should be at least 4 carbons
    carbon_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with phosphoethanolamine group and single acyl substituent at sn-1"