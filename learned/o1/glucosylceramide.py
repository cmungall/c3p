"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: CHEBI:37547 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a ceramide (sphingoid base linked to fatty acid) with 
    a glucose moiety attached via a beta-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide backbone - sphingoid base
    # Sphingoid base: Long chain amino alcohol with (E)-double bond at C4-C5 and hydroxyl groups at C1 and C3
    sphingoid_pattern = Chem.MolFromSmarts(
        "[#6][#6][#6][CH]=[CH][CH](O)[CH](N)[CH2]O"
    )
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Check for fatty acid linked via amide bond to sphingoid base
    amide_pattern = Chem.MolFromSmarts(
        "C(=O)N[CH]"
    )
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No fatty acid amide linkage found"

    # Check for glucose head group attached via beta-glycosidic bond
    # Glucose pattern (beta-D-glucose)
    glucose_pattern = Chem.MolFromSmarts(
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1"
    )
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No glucose moiety found"

    # Check for glycosidic linkage between glucose and sphingoid base
    # Glycosidic linkage between the anomeric carbon of glucose and primary hydroxyl of sphingoid base
    glycosidic_linkage_pattern = Chem.MolFromSmarts(
        "[C@@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)-[O]-[CH2]-[CH](N)-C(=O)"
    )
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No beta-glycosidic linkage to ceramide found"

    # Additional checks (optional)
    # Verify that the glucose is connected via beta linkage
    # Verify chain lengths of fatty acid and sphingoid base (usually long chains)

    return True, "Contains ceramide backbone with beta-D-glucose attached via beta-glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37547',
                              'name': 'glucosylceramide',
                              'definition': 'Any of the cerebrosides in which the monosaccharide head group is glucose.',
                              'parents': ['CHEBI:33563']},
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
        'stdout': None
    }