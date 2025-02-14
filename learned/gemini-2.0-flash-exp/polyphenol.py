"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol has 2 or more benzene rings, each with at least one directly attached -OH group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polyphenol, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    
    # Get all the benzene ring matches
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    # Count of benzene rings with at least one OH group
    oh_count = 0
    
    # Loop through all benzene rings
    for match in benzene_matches:
      benzene_atoms = [mol.GetAtomWithIdx(i) for i in match]
      has_oh = False
      # check each atom for neighbor oxygen with one hydrogen atom
      for atom in benzene_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Ensure it's an OH group (single bond to O, O has one H)
                if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                  
                  h_count = 0
                  for n2 in neighbor.GetNeighbors():
                    if n2.GetAtomicNum() == 1:
                      h_count += 1
                  if h_count == 1:
                    has_oh = True
                    break
        if has_oh:
          break
      if has_oh:
         oh_count += 1
    
    if oh_count >= 2:
       return True, "Contains two or more benzene rings, each with at least one directly attached -OH group."
    else:
      return False, f"Found {oh_count} benzene rings with directly attached -OH groups, needs at least two."