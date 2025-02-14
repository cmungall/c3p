"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a compound consisting of a hydroxysteroid backbone with sugar moieties attached via glycosidic bonds.

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
    
    # Define steroid nucleus pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('''
        [#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]([#6][#6]2)[#6][#6][#6]4[C,H]([#6][#6]3)[#6][#6][#6]([#6]4)[#6]
    ''')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for hydroxyl groups attached to the steroid backbone
    hydroxyl_pattern = Chem.MolFromSmarts('[#6][$([#6]O)]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found on the steroid backbone"
    
    # Define sugar (glycoside) pattern (pyranose ring)
    sugar_pattern = Chem.MolFromSmarts('''
        [C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[O,C]-[C,H]1
    ''')
    # Check for sugar moieties attached via glycosidic bonds
    glycosidic_bond_pattern = Chem.MolFromSmarts('[C;!R][O][C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H]1O')
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No sugar moieties attached via glycosidic bonds found"
    
    # Additional check: ensure that the sugar is attached to the steroid backbone
    # Identify connection between steroid backbone and sugar moiety
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    connected = False
    for s_match in steroid_matches:
        s_atoms = set(s_match)
        for sugar_match in sugar_matches:
            sugar_atoms = set(sugar_match)
            for bond in mol.GetBonds():
                if (bond.GetBeginAtomIdx() in s_atoms and bond.GetEndAtomIdx() in sugar_atoms) or \
                   (bond.GetBeginAtomIdx() in sugar_atoms and bond.GetEndAtomIdx() in s_atoms):
                    connected = True
                    break
            if connected:
                break
        if connected:
            break
    if not connected:
        return False, "Steroid backbone not connected to sugar moiety via glycosidic bond"

    return True, "Contains a hydroxysteroid backbone with sugar moieties attached via glycosidic bonds"