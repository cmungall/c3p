"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a steroid with a hydroxyl group at the 3-beta position.
    It has a tetracyclic ring system closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core SMARTS pattern (simplified - does not enforce stereochemistry on ring fusion)
    # The core is based on the perhydrophenanthrene ring system, three 6 membered rings and one 5 membered ring.
    # The carbon on which the hydroxy group is attached is C3.
    steroid_core_pattern = Chem.MolFromSmarts("[C]12[C]([C]3[C]([C]4[C]([C]([C]1)CC[C@@H]2C)[C](CC4)CC3)C)") # This SMARTS is much simpler than the previous version
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"
    
    # Check for a hydroxyl group at the 3-position (we have to find the 3-position and then check for an OH)
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if core_match:
        c3 = core_match[0] # This is the C connected to the hydroxy group, in core_match, this is the first C defined.
        c3_atom = mol.GetAtomWithIdx(c3)
        has_oh = False
        for neighbor in c3_atom.GetNeighbors():
             if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                has_oh = True
                break
        if not has_oh:
            return False, "No hydroxyl group at the 3-position"
    else:
        return False, "No match of steroid_core_pattern"
    
    # Check for a side chain at position 17 (simplified SMARTS for a C attached to C17).
    # It can have carbon chain attached to it. C17 is the carbon attached to the 5 membered ring in the core.
    c17 = core_match[5] # C17 in the simplified core pattern
    c17_atom = mol.GetAtomWithIdx(c17)

    has_sidechain = False
    for neighbor in c17_atom.GetNeighbors():
         if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() < 3:
             has_sidechain = True
             break
    if not has_sidechain:
        return False, "No side chain at C17"
        
    return True, "Contains a steroid core, a hydroxyl group at the 3-position, and a side chain at C17."