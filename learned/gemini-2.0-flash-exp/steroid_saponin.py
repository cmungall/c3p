"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a saponin derived from a hydroxysteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for steroid core and glycosidic linkages
    # Steroid core (gonane ring system with common substitutions)
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C]([C]([C]([C]2[C]1([C]([C]([C]3[C]2([C]([C]([C]([C]3)([C]4)[C])([C]5)[C])([C]4)C)C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found."
    
    # Check for glycosidic linkages: Carbon-O-Sugar
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[CX4][OX2][C](O)")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if len(glycosidic_matches) == 0 :
      return False, "No glycosidic linkage found. Thus not a saponin"
    
    # Check for sugar presence
    sugar_pattern = Chem.MolFromSmarts("[C][C](O)[C](O)[C](O)[C](O)[C]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0 :
      return False, "No sugar moiety found. Thus not a saponin"
    
    # Check molecular weight: steroid saponins are relatively large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
      return False, "Molecular weight too low to be a steroidal saponin"

    # Check number of carbons and oxygens: steroid saponins are large and have several oxygens.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 20:
        return False, "Too few carbons to be steroidal saponin"
    if o_count < 5:
        return False, "Too few oxygens to be steroidal saponin"

    return True, "Contains steroid core with glycosidic linkages and sugar moieties (i.e. is a steroid saponin)."