"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue linked by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the sulfonate group (-OS(=O)(=O)O-) or sulfonic acid (-S(=O)(=O)OH).
    sulfonate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX1]")
    sulfonic_acid_pattern = Chem.MolFromSmarts("[S](=O)(=O)[OX1H0]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    sulfonic_acid_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    if not sulfonate_matches and not sulfonic_acid_matches:
       return False, "No sulfonate or sulfonic acid group found."

    # Get atoms of the sulfonate group
    sulfur_atoms = []
    if sulfonate_matches:
      for match in sulfonate_matches:
        for idx in match:
          atom = mol.GetAtomWithIdx(idx)
          if atom.GetAtomicNum() == 16: # Sulfur
            sulfur_atoms.append(atom)
    elif sulfonic_acid_matches:
      for match in sulfonic_acid_matches:
        for idx in match:
          atom = mol.GetAtomWithIdx(idx)
          if atom.GetAtomicNum() == 16: # Sulfur
            sulfur_atoms.append(atom)
    
    # 2. Check the carbons directly bonded to the sulfur
    connected_carbons = []
    for sulfur_atom in sulfur_atoms:
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6: # Carbon
                connected_carbons.append(neighbor)
    if not connected_carbons:
        return False, "No carbon directly bonded to sulfur."

    # 3. Check if the carbon is part of a lipid-like moiety
    is_lipid = False
    for carbon_atom in connected_carbons:
      # Check for long chains connected to the carbon
      chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
      
      # Check if this carbon is part of such a chain
      match_list = mol.GetSubstructMatches(chain_pattern)
      for match in match_list:
          if carbon_atom.GetIdx() in match:
            is_lipid = True
            break
      if is_lipid:
        break;

      #Check for fatty acid chain
      fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3](=[OX1])")
      match_list = mol.GetSubstructMatches(fatty_acid_pattern)
      for match in match_list:
          if carbon_atom.GetIdx() in match:
            is_lipid = True
            break
      if is_lipid:
        break;

      # Check for common sugars
      sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][OX2]") #simplified sugar ring
      match_list = mol.GetSubstructMatches(sugar_pattern)
      for match in match_list:
        if carbon_atom.GetIdx() in match:
            is_lipid = True
            break
      if is_lipid:
        break;

    if not is_lipid:
        return False, "Carbon connected to sulfur is not part of a lipid-like moiety"


    # 4. Check molecular weight to discard obvious non-lipids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
      return False, "Molecular weight too low for a sulfolipid"


    return True, "Contains a sulfonic acid residue linked to a lipid via a carbon-sulfur bond."