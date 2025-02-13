"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from oxidation of C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:  # Ensuring focus on robust C20 derivatives
        return False, f"Carbon count ({c_count}) not typical for icosanoids"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Insufficient number of oxygens ({o_count}), expect oxidized features in icosanoids"

    # Check for cyclopentane ring (typical in prostaglandins)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "Expected cyclopentane ring not detected in icosanoids"

    # Check for functional groups: keto (C=O) and hydroxyl (OH)
    keto_pattern = Chem.MolFromSmarts("C=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(keto_pattern):
        return False, "Lacks keto group common in icosanoids"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Lacks hydroxyl group common in icosanoids"

    # Check for multiple double bonds (characteristic of EFAs)
    double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2])
    if double_bonds < 4:  # Must accommodate polyunsaturation typical of parent acids
        return False, f"Insufficient number of double bonds ({double_bonds}) for an icosanoid"

    # Consider ester or amide linkages typical in complex icosanoids
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0,R0]")  # Ester linkage
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")       # Amide linkage
    if mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern):
        return True, "Structure contains substantial icosanoid-linked features via ester or amide linkages"

    return True, "Molecule generally fits refined structural criteria suggestive of an icosanoid"