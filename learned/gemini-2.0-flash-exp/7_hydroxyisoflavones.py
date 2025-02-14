"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the 7-hydroxyisoflavone core structure
    # Looking for c1cc(O)cc2c(c1)c(=O)oc2 OR c1cc(O)cc2c(c1)oc(=O)c2
    core_pattern1 = Chem.MolFromSmarts('c1cc(O)cc2c(c1)c(=O)oc2')
    core_pattern2 = Chem.MolFromSmarts('c1cc(O)cc2c(c1)oc(=O)c2')
    if not (mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)) :
        return False, "Molecule does not have the 7-hydroxyisoflavone core structure"

    # Look for carbonyl group attached to a phenyl or a substituted phenyl group
    phenyl_group_pattern = Chem.MolFromSmarts('c1ccccc1C(=O)')
    phenyl_substituted_group_pattern = Chem.MolFromSmarts('c1ccccc1-!@C(=O)')
    
    if not (mol.HasSubstructMatch(phenyl_group_pattern) or mol.HasSubstructMatch(phenyl_substituted_group_pattern)):
         return False, "Isoflavone core must have carbonyl group bonded to a phenyl"   

    return True, "Molecule has the 7-hydroxyisoflavone core with a hydroxyl group at the 7-position"