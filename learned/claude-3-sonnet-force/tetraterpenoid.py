"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:38169 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 backbone) which
    may be rearranged or modified by removal of atoms (typically methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least 40 carbon atoms (tetraterpene backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 40:
        return False, "Less than 40 carbon atoms (not a tetraterpene backbone)"
    
    # Check for long carbon chains and rings (typical terpenoid skeleton)
    largest_chain = max(len(chain) for chain in Chem.Molecules.EnumerateAliphaticSpiroCycles(mol))
    if largest_chain < 10:
        return False, "No long carbon chains/rings typical of terpenoids"
    
    # Count degree of unsaturation (# of rings + # of double bonds)
    ri = rdMolDescriptors.CalcNumRotatableBonds(mol)
    dou_bonds = sum(bond.GetIsAromatic() for bond in mol.GetBonds())
    rings = mol.GetRingInfo().NumRings()
    dou_eq = rings + dou_bonds
    if dou_eq > 16:
        return False, "Too many rings/double bonds for a terpenoid"
    
    # Check for >=5 quaternary carbons (typical of terpenoids)
    quat_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and sum(bond.GetIsAromatic() for bond in atom.GetBonds()) == 0 and atom.GetDegree() == 4)
    if quat_c < 5:
        return False, "Not enough quaternary carbons for a terpenoid"
    
    # Check molecular weight (typically >500 Da for a tetraterpene)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a tetraterpene"
    
    return True, "Contains a tetraterpene backbone (>=40 carbons) with typical terpenoid features"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:38169',
        'name': 'tetraterpenoid',
        'definition': 'Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups).',
        'parents': ['CHEBI:35686', 'CHEBI:24642']
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
    'num_true_positives': 116,
    'num_false_positives': 4,
    'num_true_negatives': 182401,
    'num_false_negatives': 63,
    'num_negatives': None,
    'precision': 0.9670329670329671,
    'recall': 0.6481481481481481,
    'f1': 0.7757282307232771,
    'accuracy': 0.9995554449262074
}