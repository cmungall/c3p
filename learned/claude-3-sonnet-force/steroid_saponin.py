"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:27559 steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.

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
    
    # Look for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@@]12[C@H](C[C@H](C[C@]1([H])[H])[C@]1([H])C[C@@H](C)C[C@]2([H])[C@@]1([H])C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for glycoside groups (sugar moieties)
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    glycoside_matches = mol.GetSubstructMatches(glycoside_pattern)
    if not glycoside_matches:
        return False, "No glycoside groups found"
    
    # Count rotatable bonds to verify sugar chains
    n_rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Sugar chains too short or missing"
    
    # Check molecular weight - steroid saponins typically >500 Da
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for steroid saponin"
    
    # Look for common steroid saponin cores
    core_smarts = ["OC1(CCC2C1CCC3C2CCC4C3CCC(O)C4C)C", # Spirostan core
                   "OC1(CCC2C3CCC4C5CCC6C7CCC8C9CCC(O)C(C9)C8C7CCC6C5CCC4C3CCC12C)C"] # Cucurbitan core
    for core in core_smarts:
        core_pattern = Chem.MolFromSmarts(core)
        if mol.HasSubstructMatch(core_pattern):
            return True, "Contains steroid backbone and glycoside groups"
    
    return True, "Contains steroid backbone and glycoside groups (no common core matched)"