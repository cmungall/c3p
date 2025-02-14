"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: CHEBI:85815 1-O-acylglycerophosphoethanolamine

A glycerophosphoethanolamine having an unspecified O-acyl substituent at the 1-position of the glycerol fragment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

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

    # Look for glycerol backbone with phosphoethanolamine head group
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    pe_pattern = Chem.MolFromSmarts("P(OCCN)(O)(=O)[OX2][CHX4][CHX4][OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(glycerol_pattern) or not mol.HasSubstructMatch(pe_pattern):
        return False, "No glycerophosphoethanolamine backbone found"
    
    # Look for acyl chain (long carbon chain attached to glycerol oxygen)
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[OX2][CHX4][CHX4][OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl chain found"
    
    # Count rotatable bonds to verify long acyl chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short"

    # Check molecular weight - typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20 or o_count < 6 or n_count != 1 or p_count != 1:
        return False, "Incorrect atom counts for 1-O-acylglycerophosphoethanolamine"

    return True, "Contains glycerophosphoethanolamine backbone with acyl chain attached at 1-position"