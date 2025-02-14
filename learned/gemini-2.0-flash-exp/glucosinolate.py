"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for thioglucoside core (glucose linked to sulfur)
    thioglucoside_pattern = Chem.MolFromSmarts("[C]1[O][C]([C]([C]([C]([O]1)CO)[O])[O])[S]")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
      return False, "Thioglucoside core not found"

    # 2. Check for the central C atom bonded to S, N, and a side-chain
    central_carbon_pattern = Chem.MolFromSmarts("[S][C](=[N])")
    matches = mol.GetSubstructMatches(central_carbon_pattern)
    if not matches:
        return False, "Central carbon with S and N not found"
    
    
    #3. Check for sulfonated oxime and its stereochemistry
    sulfonated_oxime_pattern = Chem.MolFromSmarts("C=N/OS(=O)(=O)[O-]")
    sulfonated_oxime_matches = mol.GetSubstructMatches(sulfonated_oxime_pattern)
    sulfonated_oxime_pattern_nonanionic = Chem.MolFromSmarts("C=N/OS(=O)(=O)O")
    sulfonated_oxime_matches_nonanionic = mol.GetSubstructMatches(sulfonated_oxime_pattern_nonanionic)


    if not sulfonated_oxime_matches and not sulfonated_oxime_matches_nonanionic:
      return False, "Sulfonated oxime group not found or incorrect stereochemistry"

    if sulfonated_oxime_matches and sulfonated_oxime_matches_nonanionic:
      return False, "Found a sulfonated oxime group but could not determine anion status"
    
    
    
    # 4. Check the connectivity and stereochemistry between S of thioglucoside to the central C atom, and the oxime to the central C atom.
    
    
    
    for match in matches:
      center_c_index = match[0]
      
      sulfur_neighbors = [atom.GetIdx() for atom in mol.GetAtomWithIdx(center_c_index).GetNeighbors() if atom.GetAtomicNum() == 16]
      nitrogen_neighbors = [atom.GetIdx() for atom in mol.GetAtomWithIdx(center_c_index).GetNeighbors() if atom.GetAtomicNum() == 7]

      if len(sulfur_neighbors) != 1 or len(nitrogen_neighbors) != 1 :
        return False, "Central Carbon is not connected to exactly one S and one N"
      
      sulfur_index = sulfur_neighbors[0]
      nitrogen_index = nitrogen_neighbors[0]
      
      sulfur_matches = mol.GetSubstructMatches(thioglucoside_pattern)
      
      if not any(sulfur_index in match for match in sulfur_matches):
        return False, "Central sulfur not part of the thioglucoside group"
    
    if sulfonated_oxime_matches:
          return True, "Molecule matches all required patterns for glucosinolate with anionic sulfate"
    else:
        return True, "Molecule matches all required patterns for glucosinolate with non-anionic sulfate"