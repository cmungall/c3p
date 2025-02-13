"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:36287 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is defined as any monosaccharide containing an alcoholic
    hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of phosphate group (OP(O)(O)=O or [P+](O)(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("[P+](O)(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Check for sugar backbone (ring with multiple hydroxy groups)
    sugar_pattern = Chem.MolFromSmarts("[OX2]r1[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone found"
    
    # Check for ester bond between phosphate and hydroxy group
    ester_pattern = Chem.MolFromSmarts("[OX2][P+](O)(O)(=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond between phosphate and sugar"
    
    # Count phosphate groups (should be 1 for monosaccharide)
    phosphate_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphate_count != 1:
        return False, f"Found {phosphate_count} phosphate groups, expected 1"
    
    # Check molecular weight (typically <500 Da for monosaccharide)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for monosaccharide"
    
    # Count carbons and oxygens (rough check for carbohydrate)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 3 or o_count < 3:
        return False, "Insufficient carbons or oxygens for carbohydrate"
    
    return True, "Contains a sugar backbone with a phosphate ester"