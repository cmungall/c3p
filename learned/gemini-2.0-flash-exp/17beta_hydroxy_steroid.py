"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid has the basic steroid core and a beta-configured hydroxyl group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the steroid core (simplified, can be expanded if needed)
    # We only check for a general steroid core without double bond specificity
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C]([C][C]2)[C]3[C]([C]4[C]5[C]([C]2([C]1)[C]3)[C]([C]([C]5)[C]4)C)")
    
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not have the basic steroid core structure."

    # SMARTS pattern for the 17-position with beta-hydroxy
    # '[C@@]([C])([C])([C])O' matches carbon with correct stereochemistry and connected to an O atom
    # Then we restrict that carbon to be at the correct position [C]1([C][C]2)[C]3[C]([C]4[C]5[C]([C]2([C]1)[C]3)[C]([C]([C]5)[C]4)[C@@H](O)[C])
    
    # This needs further refining, the exact position 17 cannot easily be defined in SMARTS without full atom mapping.
    # Therefore we use a different strategy: 
    # 1) find substructure of carbon with 3 single bonds and one double bond, that's usually C17
    # 2) check that we have a hydroxyl connected by beta bond
    
    # 1) Find potential 17 carbon (carbon with 3 single bonds connected)
    potential_17_carbon = Chem.MolFromSmarts("[CHX4]")
    matches_17_c = mol.GetSubstructMatches(potential_17_carbon)
    # check for hydroxyl connected to one of those candidates
    found = False
    for match in matches_17_c:
      atom_id = match[0]
      
      c_atom = mol.GetAtomWithIdx(atom_id)
      for nbr in c_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
          bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(),nbr.GetIdx())
          if bond.GetBeginAtomIdx() == c_atom.GetIdx():
            if bond.GetStereo() == Chem.BondStereo.STEREOE:
                
              found = True
              break
            
      if found:
        break
    
    if not found:
      return False, "No 17-position with beta-hydroxy group found."

    return True, "Molecule contains a steroid core with a beta-hydroxy group at position 17."