"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester bond at the 3-hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for steroid core
    # Using a SMARTS pattern that matches the core fused rings. The pattern
    # matches the three 6 membered rings plus one 5 membered ring, with possible
    # substitutions
    steroid_core_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1]2[CR1][CR1][CR1]3[CR1][CR1][CR1][CR1]4[CR1][CR1]3[CR1]214")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    # 2. Check for hydroxyl group at position 3
    # We need to find a carbon in the steroid core that is bonded to an oxygen, which is also the 3rd carbon of the A ring.
    # This approach assumes a standard steroid numbering
    # The pattern below matches the carbon at position 3 of the sterol
    sterol_C3_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CX4]([OX2])[CR1]2[CR1][CR1][CR1]3[CR1][CR1][CR1][CR1]4[CR1][CR1]3[CR1]214")

    if not mol.HasSubstructMatch(sterol_C3_pattern):
      return False, "No hydroxyl at C3 position"

    # 3. Check for ester bond connecting to the C3 oxygen
    # the pattern looks for an O linked to C3 that is part of an ester group
    ester_pattern = Chem.MolFromSmarts("[CX4]([OX2])[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    sterol_ester_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CX4]([OX2][CX3](=[OX1])[#6])[CR1]2[CR1][CR1][CR1]3[CR1][CR1][CR1][CR1]4[CR1][CR1]3[CR1]214")

    if not mol.HasSubstructMatch(sterol_ester_pattern):
        return False, "No ester bond at position 3"


    #4. check for presence of other carboxylic acids,
    # looking for any -C(=O)OH group that is not part of the ester bound to the sterol.

    free_acid_pattern = Chem.MolFromSmarts("C(=O)O[H]")

    acid_matches = mol.GetSubstructMatches(free_acid_pattern)
    if acid_matches:
          # get the atoms matched from pattern 3
      ester_atoms_matches = mol.GetSubstructMatches(sterol_ester_pattern)
      if ester_atoms_matches:
        # get the oxygen atoms forming the ester bond (atom with index 1 in sterol_ester_pattern)
          ester_oxygens = set([x[1] for x in ester_atoms_matches])
      
        #check if any of the carboxylic acid oxygens are also an ester oxygen and keep only the free acids
          free_acids = [match for match in acid_matches if match[1] not in ester_oxygens]
          if free_acids:
            return False, "Free carboxylic acid group present, not a simple sterol ester"
      else: #if we do not have ester at 3 position but we do have other carboxyls, then this is false.
        return False, "Free carboxylic acid group present and no ester at 3 position"


    return True, "Sterol ester detected: Contains a steroid core with an ester bond at position 3"