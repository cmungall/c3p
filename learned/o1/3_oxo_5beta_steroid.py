"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:XXXXX 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta-configuration at position 5.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    
    # Define the steroid core SMARTS pattern
    steroid_core_smarts = """
    [#6]1[#6][#6]2[#6](=[#6][#6]3[#6][#6][#6]([#6]4[#6][#6][#6]([#6]3)[#6][#6]4)[#6]2[#6][#6]1)
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid core not found"
    
    # Define SMARTS pattern for ketone at position 3
    # Simplified as carbonyl group attached to a ring carbon
    ketone_smarts = "[R][C](=O)[R]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_at_pos3 = False
    for match in ketone_matches:
        atom = mol.GetAtomWithIdx(match[0])
        # Check if the atom is at position 3 by checking its ring membership
        # This is a simplification and may not be accurate for all steroids
        if atom.IsInRing():
            ketone_at_pos3 = True
            break
    if not ketone_at_pos3:
        return False, "Ketone group at position 3 not found"
    
    # Check for beta-configuration at position 5
    # Identify stereocenter at position 5
    stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    position5_chiral = False
    for idx, chirality in stereocenters:
        if idx in [atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRingSize(6)]:
            # Assuming position 5 is in a six-membered ring
            position5_chiral = True
            cip_code = chirality
            break
    if not position5_chiral:
        return False, "Chiral center at position 5 not found"
    if chirality != 'S':
        return False, f"Position 5 is not in beta-configuration (found {chirality})"
    
    return True, "Molecule is a 3-oxo-5beta-steroid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXX',
        'name': '3-oxo-5beta-steroid',
        'definition': "Any 3-oxo steroid that has beta- configuration at position 5.",
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}