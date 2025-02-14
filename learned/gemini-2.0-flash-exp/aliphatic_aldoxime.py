"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is an aldoxime where the carbon attached to the oxime group is derived from an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for the oxime group (C=N-O)
    oxime_pattern = Chem.MolFromSmarts("[CH1]=[N]-O")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No aldoxime group found"
    
    # Get the carbon atoms from the oxime group
    matches = mol.GetSubstructMatches(oxime_pattern)

    # Function to check if the chain connected to the oxime carbon is aliphatic
    def is_aliphatic_chain(atom, visited_atoms):
      visited_atoms.add(atom.GetIdx())
      for neighbor in atom.GetNeighbors():
          if neighbor.GetIdx() not in visited_atoms: # Check for visited atoms to prevent loops
            if neighbor.GetIsAromatic():
                return False
            #Recursive step
            if neighbor.GetAtomicNum() == 6:
                if not is_aliphatic_chain(neighbor, visited_atoms):
                    return False
      return True

    # 3. Check for aliphatic nature of the chain connected to the carbon
    for match in matches:
        carbon_index = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        if not is_aliphatic_chain(carbon_atom,set()):
            return False, "The carbon of the oxime group is not part of an aliphatic chain"

    # 4. Return True if all conditions are met
    return True, "Molecule is an aliphatic aldoxime"