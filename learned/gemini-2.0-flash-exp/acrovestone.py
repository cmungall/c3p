"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone based on its SMILES string.
    Acrovestones are polyphenols, often with an isoflavone or flavone backbone,
    multiple hydroxyl groups, and often glycosylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for isoflavone or flavone core
    isoflavone_pattern = Chem.MolFromSmarts("c1cc(cc(c1)C2=CC(=O)c3ccccc3O2)O") # basic isoflavone or flavone core
    if not mol.HasSubstructMatch(isoflavone_pattern):
       isoflavone_pattern = Chem.MolFromSmarts("c1cc(cc(c1)C2=CC(=O)c3ccc(O)cc3O2)O") # basic isoflavone or flavone core with two O
       if not mol.HasSubstructMatch(isoflavone_pattern):
         return False, "No isoflavone/flavone core found"

    # 2. Check for multiple hydroxyl groups (at least 2, often more)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    if hydroxyl_count < 2:
        return False, f"Too few hydroxyl groups: {hydroxyl_count}"
    
    # 3. Check for glycosylation (at least 1 sugar moiety)
    sugar_pattern = Chem.MolFromSmarts("[C](-[O])-[C](-[O])-[C](-[O])-[C](-[O])-[C](-[O])-[C](-[O])") # check for 6-member sugar
    
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    
    
    sugar_pattern_pentose = Chem.MolFromSmarts("[C](-[O])-[C](-[O])-[C](-[O])-[C](-[O])-[C](-[O])") #check for 5-member sugar
    pentose_matches = mol.GetSubstructMatches(sugar_pattern_pentose)


    if len(sugar_matches) == 0 and len(pentose_matches) == 0:
      
      is_glycoside = False
    else:
      is_glycoside = True

    if not is_glycoside and hydroxyl_count<3 : #some acrovestones are not glycosides but tend to have more OHs
          return False, "No glycosylation detected, must have 3 or more hydroxyl groups"


    # Check for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("[OCH3]")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))

    #4 check for prenyl pattern

    prenyl_pattern = Chem.MolFromSmarts("CC(C)(O)CCC")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)

    # # 5. Verify substitution at expected positions (7, 4', 5 or 6) - This is difficult to do in a generalized way, not checking for now

    #additional checks
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Too few carbons for acrovestone"

    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
         return False, "Too few oxygens for acrovestone"
    
    return True, "Matches acrovestone structural features"