"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:26216 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units of uronic acids and glycosamines,
    often partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for uronic acid and glycosamine patterns
    uronic_acid_pattern = Chem.MolFromSmarts("[C&R1](=O)[O-]")
    glycosamine_pattern = Chem.MolFromSmarts("[N;R1][C;R2][C;R2][C;R2][C;R2][O;R1]")
    
    uronic_acids = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamines = mol.GetSubstructMatches(glycosamine_pattern)
    
    # Check for alternating units
    if not uronic_acids or not glycosamines:
        return False, "No uronic acid or glycosamine units found"
    
    # Count sulfates
    sulfate_pattern = Chem.MolFromSmarts("[S+4](=O)(=O)(-[O-])")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    n_sulfates = len(sulfate_matches)
    
    # Check for polysaccharide structure
    n_atoms = mol.GetNumAtoms()
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_atoms < 20 or n_rotatable < 5:
        return False, "Molecule too small to be a polysaccharide"
    
    # Classify based on structure
    if n_sulfates > 0:
        return True, f"Contains uronic acid and glycosamine units, {n_sulfates} sulfate groups"
    else:
        return True, "Contains uronic acid and glycosamine units in alternating pattern"