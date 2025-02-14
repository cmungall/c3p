"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the chromene core with flexibility for saturation, and a connection point.
    # [CH2X4,CHX3] means either CH2 (saturated ring) or CH (unsaturated)
    # [*] means that this atom must have an extra attachment to something.
    chromene_core_pattern = Chem.MolFromSmarts('c1ccccc1[CH2X4,CHX3]([*])O[CHX4,CHX3]')
    if not mol.HasSubstructMatch(chromene_core_pattern):
          return False, "Does not contain a 1-benzopyran core"

    # Define aryl substituent, connected via a single bond to the core. 
    #  Allow for substituted aryl rings at position 4.
    aryl_substituent_pattern = Chem.MolFromSmarts('[c]1:[c]:[c]:[c]:[c]:[c]1')
    
    # Find the attachment point to the ring.
    matches = mol.GetSubstructMatches(chromene_core_pattern)
    
    found_aryl_at_pos_4 = False
    for match in matches:
      core_C4_index = match[1] # the position bonded to the aryl.
      # iterate over the bonds of the matching atom.
      for bond in mol.GetAtomWithIdx(core_C4_index).GetBonds():
        if bond.GetOtherAtomIdx(core_C4_index) == match[2]:
          # skip the core bond
          continue
        other_atom = mol.GetAtomWithIdx(bond.GetOtherAtomIdx(core_C4_index))
        if mol.HasSubstructMatch(aryl_substituent_pattern,[other_atom.GetIdx()]):
          found_aryl_at_pos_4 = True
          break
      if found_aryl_at_pos_4:
        break
          
    if not found_aryl_at_pos_4:
        return False, "Does not have an aryl substituent at position 4 of the 1-benzopyran core."

    return True, "Contains a 1-benzopyran core with an aryl substituent at position 4."