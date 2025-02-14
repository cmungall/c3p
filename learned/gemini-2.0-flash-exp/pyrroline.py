"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:35419 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a 5-membered ring with one nitrogen and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Patterns to match pyrroline core with the double bond in different locations and the nitrogen having sp2 or sp3 character.
    # Pattern 1 : N has 1 double bond within the ring and 1 single bond
    pyrroline_pattern1 = Chem.MolFromSmarts("[NX2]1-[#6]=[#6]-[#6]-[#6]1")
    # Pattern 2 : N has 2 single bonds
    pyrroline_pattern2 = Chem.MolFromSmarts("[NX3]1-[#6]=[#6]-[#6]-[#6]1")
    pyrroline_pattern3 = Chem.MolFromSmarts("[NX3]1-[#6]-[#6]=[#6]-[#6]1")
    pyrroline_pattern4 = Chem.MolFromSmarts("[NX3]1-[#6]-[#6]-[#6]=[#6]1")


    if not (mol.HasSubstructMatch(pyrroline_pattern1) or mol.HasSubstructMatch(pyrroline_pattern2) or mol.HasSubstructMatch(pyrroline_pattern3) or mol.HasSubstructMatch(pyrroline_pattern4)):
       return False, "No pyrroline core structure found"
    
    #Ensure there is only one double bond within the ring and it is part of the pyrroline core
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    matches = mol.GetSubstructMatches(double_bond_pattern)
    ring_matches = 0
    for match in matches:
      is_ring = mol.GetAtomWithIdx(match[0]).IsInRing() and mol.GetAtomWithIdx(match[1]).IsInRing()
      if is_ring:
        #Check if the double bond is part of the pyrroline core.
        is_pyrroline_bond = False
        substructure_matches = mol.GetSubstructMatches(pyrroline_pattern1)
        for substructure_match in substructure_matches:
          if match[0] in substructure_match and match[1] in substructure_match:
              is_pyrroline_bond = True
              break
        if is_pyrroline_bond == False:
          substructure_matches = mol.GetSubstructMatches(pyrroline_pattern2)
          for substructure_match in substructure_matches:
             if match[0] in substructure_match and match[1] in substructure_match:
               is_pyrroline_bond = True
               break
        if is_pyrroline_bond == False:
          substructure_matches = mol.GetSubstructMatches(pyrroline_pattern3)
          for substructure_match in substructure_matches:
            if match[0] in substructure_match and match[1] in substructure_match:
              is_pyrroline_bond = True
              break
        if is_pyrroline_bond == False:
            substructure_matches = mol.GetSubstructMatches(pyrroline_pattern4)
            for substructure_match in substructure_matches:
                if match[0] in substructure_match and match[1] in substructure_match:
                  is_pyrroline_bond = True
                  break
        
        if is_pyrroline_bond:
            ring_matches += 1

    if ring_matches != 1:
      return False, "The pyrroline ring must have exactly 1 double bond."

    # Check that this is not a pyrrole, since a pyrrole contains 2 double bonds within the ring.
    pyrrole_pattern = Chem.MolFromSmarts("[nX2]1[cX3]=[cX2][cX2]=[cX2]1")
    if mol.HasSubstructMatch(pyrrole_pattern):
        return False, "It's a pyrrole not a pyrroline"

    return True, "Pyrroline structure detected"