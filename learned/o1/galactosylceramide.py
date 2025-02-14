"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside where the monosaccharide head group is galactose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the galactose moiety SMARTS pattern
    galactose_smarts = '[C@H]1(O)[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1'  # β-D-galactose
    galactose_pattern = Chem.MolFromSmarts(galactose_smarts)
    if galactose_pattern is None:
        return False, "Could not parse galactose SMARTS pattern"

    # Check for galactose moiety
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose moiety found"

    # Define the ceramide backbone SMARTS pattern
    ceramide_smarts = 'O=C(N[C@@H](CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)CCCCC)CCCCCCCCCCCCCCCCCC'  # Simplified ceramide structure
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if ceramide_pattern is None:
        return False, "Could not parse ceramide SMARTS pattern"

    # Check for ceramide backbone
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for amide bond between fatty acid and sphingosine base
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if amide_pattern is None:
        return False, "Could not parse amide bond SMARTS pattern"

    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Check for long hydrocarbon chains
    aliphatic_chain_pattern = Chem.MolFromSmarts('CCCCCCCCCCCCCCCCCCCC')  # 20 carbons
    if aliphatic_chain_pattern is None:
        return False, "Could not parse aliphatic chain SMARTS pattern"

    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "No long aliphatic chains found"

    # Verify stereochemistry of sphingosine base
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 3:
        return False, f"Expected at least 3 chiral centers, found {len(chiral_centers)}"

    # If all checks pass
    return True, "Molecule is classified as a galactosylceramide"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'galactosylceramide',
        'definition': 'Any of the cerebrosides in which the monosaccharide head group is galactose.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}