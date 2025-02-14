"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:28053 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and an aldehyde group at one end,
    which can exist in open-chain or cyclic forms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate exact molecular weight for validation
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 180:
        return False, f"Molecular weight {mol_wt:.2f} is not within typical range for aldopentoses"

    # Count number of carbon atoms
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(c_atoms)
    if num_carbons != 5:
        return False, f"Contains {num_carbons} carbon atoms, should be 5 for a pentose"

    # Check for aldehyde group in open-chain form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H][#6;X3](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains aldehyde group indicative of aldopentose in open-chain form"

    # Check for hemiacetal ring (cyclic form)
    hemiacetal_pattern = Chem.MolFromSmarts("[C@H]1([O])[O][C@@H]([O])[C@H]([O])[C@@H]1[O]")
    if mol.HasSubstructMatch(hemiacetal_pattern):
        return True, "Contains cyclic hemiacetal form indicative of aldopentose"

    # Check for furanose ring (5-membered ring with oxygen)
    furanose_pattern = Chem.MolFromSmarts("C1OC([OH])C([OH])C1[OH]")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains furanose ring indicative of aldopentose"

    # Check for pyranose ring (6-membered ring with oxygen)
    pyranose_pattern = Chem.MolFromSmarts("C1OC([OH])C([OH])C([OH])C1[OH]")
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains pyranose ring indicative of aldopentose"

    # If none of the patterns match
    return False, "Does not match patterns for aldopentose structures"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28053',
                              'name': 'aldopentose',
                              'definition': 'A pentose with a (potential) aldehyde group at one end.',
                              'parents': ['CHEBI:16842', 'CHEBI:24869']},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}