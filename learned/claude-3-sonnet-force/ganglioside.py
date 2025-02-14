"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: CHEBI:17619 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid with one or more sialic acids linked on the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ceramide backbone (long alkyl chain + amide)
    ceramide_pattern = Chem.MolFromSmarts("CCCCC(=O)NC")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"
    
    # Look for sialic acid residue(s)
    sialic_acid_pattern = Chem.MolFromSmarts("C1C(C(C(O1)CO)O)OC")
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if not sialic_acid_matches:
        return False, "No sialic acid residues found"
    
    # Look for glycosidic linkages between sugar residues
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H][C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkages found"
    
    # Count heavy atoms - gangliosides are typically large molecules
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 50:
        return False, "Too few heavy atoms for ganglioside"
    
    return True, "Contains ceramide backbone with sialic acid residue(s) and glycosidic linkages"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17619',
        'name': 'ganglioside',
        'definition': 'A glycosphingolipid containing one or more sialic acid residues.',
        'parents': ['CHEBI:18060', 'CHEBI:36374']
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 0,
    'num_true_negatives': 182438,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9787234042553192,
    'f1': 0.9892473118279569,
    'accuracy': 0.9999947166863722
}