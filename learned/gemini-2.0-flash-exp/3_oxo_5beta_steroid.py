"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    
    A 3-oxo-5beta-steroid has:
        - a steroid core
        - a ketone at the 3 position
        - a beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Steroid core (flexible substitution)
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C]([C])[C]3[C]2[C]([C])[C][C]4[C]3[C]([C])[C][C][C]4")

    # 3-oxo group (ketone at position 3)
    oxo_group_pattern = Chem.MolFromSmarts("[C]1[C](=O)[C][C]2[C]1[C]")

    # 5beta configuration at position 5 (implicit H)
    # Since position 5 is part of the steroid core, need to combine the steroid pattern with the 5beta bond.
    # The 5th atom from steroid_core_pattern is the 5th atom in a steroid A ring
    # I will use the index of the matching atoms to look for a C@@H pattern
    
    # Check for steroid core
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Molecule does not contain a steroid core"
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "Molecule does not contain a ketone at position 3"
    
    # Check for 5beta configuration
    # Get matches from the steroid core pattern
    matches = mol.GetSubstructMatches(steroid_core_pattern)
    
    # Ensure we have at least one match
    if len(matches) == 0:
      return False, "No steroid core match found (should not happen)"

    # Get the atoms in the steroid core. Get the 5th one.
    match = matches[0]
    atom_5 = mol.GetAtomWithIdx(match[4])
    
    # Get neighbors of atom 5
    neighbors = atom_5.GetNeighbors()
    
    # Iterate through the neighbors
    has_beta_hydrogen = False
    for neighbor in neighbors:
        if neighbor.GetSymbol() == "H":
            # Check if bond is 'up' (implicit stereochemistry)
            bond = mol.GetBondBetweenAtoms(atom_5.GetIdx(), neighbor.GetIdx())
            if bond.GetBondDir() == Chem.BondDir.BEGINWEDGE:
                has_beta_hydrogen = True
                break
            elif bond.GetBondDir() == Chem.BondDir.NONE:
                # Look at the atom itself and get its chiral tag.
                # If the chiral tag is C@@ and it has an implicit Hydrogen, that is what we are looking for.
                if atom_5.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and atom_5.GetTotalNumHs() == 1:
                  has_beta_hydrogen = True
                  break


    if not has_beta_hydrogen:
      return False, "Hydrogen at position 5 is not in beta configuration"

    return True, "Molecule is a 3-oxo-5beta-steroid"