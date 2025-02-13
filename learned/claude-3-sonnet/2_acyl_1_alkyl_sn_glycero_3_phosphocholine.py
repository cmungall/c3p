"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for complete phosphocholine group structure
    phosphocholine = Chem.MolFromSmarts("[O][P](=[O])([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "No phosphocholine group found"
    
    # Count phosphocholine groups - should be exactly one
    phosphocholine_matches = len(mol.GetSubstructMatches(phosphocholine))
    if phosphocholine_matches != 1:
        return False, f"Found {phosphocholine_matches} phosphocholine groups, need exactly 1"

    # Check for glycerol backbone with correct stereochemistry and substitution pattern
    # [C@H] ensures the sn-glycero configuration
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[#6].[CH2X4]-[OX2]-[PX4].[C@H](-[OX2]-[CX3]=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing correct glycerol backbone structure with required substitution pattern"

    # Check for ether linkage at position 1 (any alkyl)
    ether_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[CX4]")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No alkyl ether found at position 1"

    # Check for acyl group at position 2 (includes formyl)
    acyl_pattern = Chem.MolFromSmarts("[OX2]-[CX3](=[OX1])")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl group found at position 2"

    # Verify correct connectivity
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "Incorrect connectivity in glycerol backbone"

    # Check stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol)
    if not any('@' in smiles for match in chiral_centers):
        return False, "Missing or incorrect stereochemistry"

    # Verify minimum size requirements
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 8:  # Minimum size considering the basic structure
        return False, "Molecule too small for expected structure"

    # Check for required number of oxygens
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 6:  # 1 ether + 2 acyl + 3 phosphate
        return False, "Insufficient oxygen atoms for required functional groups"

    return True, "Contains glycerol backbone with alkyl ether at position 1, acyl group at position 2, and phosphocholine at position 3"