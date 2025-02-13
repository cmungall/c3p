"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphocholine group (-PO4-CH2-CH2-N+(CH3)3)
    phosphocholine_pattern = Chem.MolFromSmarts("[P](=[O])([O-])(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Look for glycerol backbone with specific stereochemistry
    # [C@H] or [C@@H] indicates the specific stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H,C@@H][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry found"

    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Verify negative charge on phosphate
    negative_charge_pattern = Chem.MolFromSmarts("[O-]P(=O)")
    if not mol.HasSubstructMatch(negative_charge_pattern):
        return False, "No negative charge on phosphate group"

    # Verify positive charge on choline
    positive_charge_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")
    if not mol.HasSubstructMatch(positive_charge_pattern):
        return False, "No positive charge on choline group"

    # Count carbons and oxygens to ensure reasonable molecule size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:  # Minimum size considering glycerol, phosphocholine and small acyl chains
        return False, "Too few carbons for a phosphatidylcholine"
    if o_count < 8:  # 2 esters (4 O), 1 phosphate (4 O)
        return False, "Too few oxygens for a phosphatidylcholine"

    # Check that fatty acid chains are attached to glycerol backbone
    # This is a simplified check for carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    return True, "Contains glycerol backbone with 2 fatty acid chains and phosphocholine group with correct charges"