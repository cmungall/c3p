"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:27195 amino sugar
Any sugar having one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino group(s)
    amino_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() < 2)
    if amino_groups == 0:
        return False, "No amino group found"
    
    # Check for sugar backbone pattern
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][OX2]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone pattern found"
    
    # Check for amino group(s) attached to the sugar backbone
    amino_sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][NX3][CX4][CX4][CX4][OX2]")
    if not mol.HasSubstructMatch(amino_sugar_pattern):
        return False, "Amino group not attached to sugar backbone"
    
    # Count hydroxy and amino groups
    hydroxy_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxy_groups + amino_groups < 3:
        return False, "Too few hydroxy and amino groups for amino sugar"
    
    return True, "Molecule contains a sugar backbone with at least one hydroxy group replaced by an amino group"