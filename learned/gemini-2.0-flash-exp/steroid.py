"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is a molecule with a cyclopenta[a]phenanthrene skeleton,
    typically with methyl groups at C-10 and C-13.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the basic steroid skeleton using SMARTS
    # This SMARTS pattern tries to capture the four fused rings with flexible bonds
    # [C]12[C]([C]3[C]([C]4[C]([C]([C]1)([C]2)[C]5)([C]3)[C]4)[C]5)  
    # However, this rigid representation causes a large number of false negatives, so we focus on the carbon atoms of the ring system only
    steroid_core_smarts = "[C]1~[C]~[C]~[C]2~[C]~[C]3~[C]~[C]4~[C]~[C]~[C]~[C](~[C]1)~[C]~[C]~4~2~3"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    if not mol.HasSubstructMatch(steroid_core_pattern):
      return False, "No steroid tetracyclic ring system found"

    # Check for methyl groups at C-10 and C-13. 
    # We first need to find the numbering of atoms according to the steroid skeleton
    match = mol.GetSubstructMatch(steroid_core_pattern)

    if not match:
        return False, "Failed to get substructure match for numbering."

    # Define substructure for C10 and C13 methyls
    methyl_pattern = Chem.MolFromSmarts("[C]-C")

    # Check number of methyl groups
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)

    if len(methyl_matches) < 2:
        return False, f"Found {len(methyl_matches)} methyl groups, at least two are required"

    # Get the atoms in the ring system
    atoms = [mol.GetAtomWithIdx(x) for x in match]
    
    # Get the neighbors of each ring atom that are not in the ring system
    neighbors = []
    for i,atom in enumerate(atoms):
      for neighbor in atom.GetNeighbors():
         if neighbor.GetIdx() not in match:
            neighbors.append((i,neighbor))

    methyl_atoms = []
    for (i,atom) in neighbors:
      if atom.GetSymbol() == 'C' and len(atom.GetNeighbors())==1:
        methyl_atoms.append( (i,atom) )

    # We check for methyl groups only, assuming two of them will be at C-10 and C-13 in most cases
    if len(methyl_atoms) < 2:
        return False, f"Found {len(methyl_atoms)} methyl groups on the ring system, at least two are required."

    return True, "Contains the steroid tetracyclic ring system with at least two methyl groups"