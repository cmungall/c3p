"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: CHEBI:17855 phosphatidylethanolamine
A class of glycerophospholipids in which a phosphatidyl group is esterified to the hydroxy group of ethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for ethanolamine group
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"
    
    # Exclude molecules with choline or similar groups
    choline_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if mol.HasSubstructMatch(choline_pattern):
        return False, "Contains choline or similar group"
    
    # Look for at least two ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    elif len(ester_matches) > 2:
        return False, f"Found {len(ester_matches)} ester groups, need at most 2"
    
    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acid chains too short"
    
    # Count carbons, oxygens, nitrogens, and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidylethanolamine"
    if o_count < 6:
        return False, "Too few oxygens for phosphatidylethanolamine"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus"
    
    return True, "Contains glycerol backbone with phosphate group esterified to ethanolamine and two fatty acid chains"