"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase and a ribose or deoxyribose sugar linked by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise.
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ribose/deoxyribose sugar pattern (5-membered ring with 1 oxygen, 4 carbons, and at least 2 OH groups )
    # Consider both ribose and deoxyribose
    sugar_pattern = Chem.MolFromSmarts("[CX4][CX4]([OX2])[CX4][CX4][CX4]") # five membered ring with one O
    
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose or deoxyribose sugar found."

    # Nucleobase patterns (purine or pyrimidine)
    purine_pattern = Chem.MolFromSmarts("n1cnc2[nH]c(=[NX3])nc2c1")
    pyrimidine_pattern = Chem.MolFromSmarts("n1cc(=[OX1])[nH]c(=[NX3])c1")
    
    purine_matches = mol.GetSubstructMatches(purine_pattern)
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)

    if not purine_matches and not pyrimidine_matches:
        return False, "No purine or pyrimidine nucleobase found."

    # Check for glycosidic bond (C-N bond linking sugar and base)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[NX3]") 
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) <1:
        return False, "No C-N glycosidic bond found between sugar and base"

    # Check that the sugar and base are actually connected to each other
    sugar_atoms = set()
    for match in sugar_matches:
      for a_idx in match:
        sugar_atoms.add(a_idx)
    
    base_atoms = set()
    if purine_matches:
      for match in purine_matches:
        for a_idx in match:
          base_atoms.add(a_idx)
    if pyrimidine_matches:
        for match in pyrimidine_matches:
          for a_idx in match:
            base_atoms.add(a_idx)

    glycosidic_atoms = set()
    for match in glycosidic_matches:
        for a_idx in match:
          glycosidic_atoms.add(a_idx)

    connected = False
    for a_idx in glycosidic_atoms:
      neighbor_atoms = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(a_idx).GetNeighbors()]
      for neighbor in neighbor_atoms:
          if neighbor in sugar_atoms and neighbor in base_atoms:
              connected = True
              break
      if connected:
          break
    if not connected:
        return False, "Sugar and base not connected via glycosidic bond"

    return True, "Contains a nucleobase linked to a ribose or deoxyribose sugar."