"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17937 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is a carbohydrate that is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acyclic structure
    if AllChem.IsCyclicMolecule(mol):
        return False, "Molecule is cyclic, alditols are acyclic"
    
    # Look for repeated pattern HOCH2[CH(OH)]nCH2OH
    alditol_pattern = Chem.MolFromSmarts("[CH2O][CH](O)[CH2](O)[CH2](O)[CH2](O)[CH2](O)")
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "Missing alditol backbone pattern"
    
    # Count hydroxy groups
    num_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetHybridization() == Chem.HybridizationType.SP3)
    if num_hydroxy < 4:
        return False, "Too few hydroxy groups for alditol"
    
    # Count carbon atoms
    num_carbon = mol.GetNumAtoms(onlyObsMask=Chem.AtomProp.Query(["#6"])) 
    if num_carbon < 4 or num_carbon > 10:
        return False, "Carbon atom count outside typical range for alditols"

    # Check for absence of other heteroatoms
    het_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "O"])
    if het_atom_count > 0:
        return False, "Molecule contains heteroatoms other than C, H, O"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 76 or mol_wt > 200:
        return False, "Molecular weight outside typical range for alditols"

    return True, "Acyclic polyol matching the formula HOCH2[CH(OH)]nCH2OH"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17937',
        'name': 'alditol',
        'definition': 'A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH (formally derivable from an aldose by reduction of the carbonyl group).',
        'parents': ['CHEBI:36973', 'CHEBI:16646']
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
    'num_true_positives': 288,
    'num_false_positives': 5,
    'num_true_negatives': 182449,
    'num_false_negatives': 13,
    'num_negatives': None,
    'precision': 0.9831460674157304,
    'recall': 0.9567567567567568,
    'f1': 0.9697694179006906,
    'accuracy': 0.9998524185023272
}