"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a polyprenol chain esterified with a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str). True if molecule is a polyprenol phosphate, False otherwise, along with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for phosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("O[P](=O)(O)O"),  # monophosphate
        Chem.MolFromSmarts("O[P](=O)(O)O[P](=O)(O)O"), #diphosphate
        Chem.MolFromSmarts("O[P](=O)(O)O[P](=O)(O)O[P](=O)(O)O")  # triphosphate
    ]
    
    #Check for the presense of at least one of those patterns.
    found_phosphate = False
    phosphate_match = None
    for pattern in phosphate_patterns:
      if mol.HasSubstructMatch(pattern):
         found_phosphate = True
         phosphate_match = pattern
         break

    if not found_phosphate:
         return False, "No phosphate group found"

    #Find the polyprenol chain via the isoprenoid pattern
    isoprenoid_pattern = Chem.MolFromSmarts("CC(=C)CC")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)
    if len(isoprenoid_matches) < 3:
          return False, "Not enough isoprenoid units in chain"
    
    #Check for connectivity between the phosphate group and the polyprenol
    #We look for the substructure C-O-P
    connectivity_pattern = Chem.MolFromSmarts("[CX4,CX3]O[P]")
    connectivity_matches = mol.GetSubstructMatches(connectivity_pattern)
    
    if len(connectivity_matches) < 1:
         return False, "Phosphate group is not connected to the chain via an oxygen"
         
    # Check the length of the polyprenol chain by looking at the number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
      return False, "Chain too short for polyprenol"

    return True, "Contains a polyprenol chain esterified with a phosphate group"