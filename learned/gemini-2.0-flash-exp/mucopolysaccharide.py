"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely a mucopolysaccharide based on its SMILES string using heuristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if likely a mucopolysaccharide, False otherwise, along with a reason.
               Returns (None, None) if SMILES processing fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None # Indicate processing failure

    # Look for carboxyl groups (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    #Look for sulfate groups (-OS(=O)(=O)O or -S(=O)(=O)O)
    sulfate_pattern1 = Chem.MolFromSmarts("OS(=O)(=O)O")
    sulfate_pattern2 = Chem.MolFromSmarts("S(=O)(=O)O")
    sulfate_matches1 = mol.GetSubstructMatches(sulfate_pattern1)
    sulfate_matches2 = mol.GetSubstructMatches(sulfate_pattern2)
    sulfate_matches = sulfate_matches1 + sulfate_matches2

    # Look for amino groups and sugars 
    amino_pattern = Chem.MolFromSmarts("[NX3]")
    oxygen_pattern = Chem.MolFromSmarts("[OX2]") #only non-carbonyl oxygens
    carbon_pattern = Chem.MolFromSmarts("[CX4]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    carbon_matches = mol.GetSubstructMatches(carbon_pattern)
    
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10:
        return False, "Too small for mucopolysaccharide"

    if len(carboxyl_matches) == 0 and len(sulfate_matches) == 0:
      return False, "No carboxyl or sulfate groups, not likely a mucopolysaccharide"
    
    #Check for minimum number of O atoms and carbons
    if len(oxygen_matches) < 5:
        return False, "Too few oxygen atoms"
    if len(carbon_matches) < 10:
        return False, "Too few carbons"
    
    #Check the presence of any amino groups: these could be aminosugars, or just a side chain
    if len(amino_matches) == 0:
        return False, "No amino groups"
    
    return True, "Contains carboxyl and/or sulfate groups, amino groups, and multiple oxygen and carbon atoms"