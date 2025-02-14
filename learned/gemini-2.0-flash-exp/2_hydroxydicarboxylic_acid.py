"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxyl group on the alpha carbon
    relative to a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True and a reason string if the molecule is a 2-hydroxydicarboxylic acid,
                         False and a reason string otherwise, or None, None if input is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # 1. Check for two carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 2:
       return False, f"Molecule does not have 2 carboxyl groups. Found {len(carboxyl_matches)}"

    # 2. Check for a hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "Molecule does not contain any hydroxyl group."

    # 3. Check for at least one hydroxyl on an alpha carbon relative to *any* carboxyl
    alpha_hydroxyl_found = False
    for carboxyl_match in carboxyl_matches:
      carboxyl_carbon_idx = carboxyl_match[0] # Get index of the carboxyl carbon

      carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

      for neighbor in carboxyl_carbon.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            # Check if neighbor is a carbon
            if neighbor.GetAtomicNum() != 6:
                continue

            # Check if any of the hydroxyl groups are attached to the neighbor carbon
            for hydroxyl_match in hydroxyl_matches:
                  hydroxyl_oxygen_idx = hydroxyl_match[0]
                  hydroxyl_oxygen = mol.GetAtomWithIdx(hydroxyl_oxygen_idx)

                  hydroxyl_carbon_neighbors = [atom.GetIdx() for atom in hydroxyl_oxygen.GetNeighbors()]
                  if neighbor_idx in hydroxyl_carbon_neighbors:
                        alpha_hydroxyl_found = True
                        break
            if alpha_hydroxyl_found:
                  break

      if alpha_hydroxyl_found:
          break
    
    if not alpha_hydroxyl_found:
      return False, "Hydroxyl group is not on the alpha carbon relative to any carboxylic group."


    #4 Success
    return True, "Molecule is a 2-hydroxydicarboxylic acid"