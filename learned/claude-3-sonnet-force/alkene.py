"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:16236 alkene
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond
    and the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for double bonds
    double_bonds = mol.GetBonds(Chem.rdchem.BondType.DOUBLE)
    num_double_bonds = len(double_bonds)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds, alkenes must have exactly 1"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Alkenes must be acyclic (no rings)"

    # Check for elements other than C and H
    non_chatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    if non_chatoms:
        non_ch_symbols = [Chem.PeriodicTable.GetElementSymbol(atom.GetAtomicNum()) for atom in non_chatos]
        return False, f"Alkenes can only contain C and H, found {', '.join(non_ch_symbols)}"

    # Check molecular formula
    formula = AllChem.CalcMolFormula(mol)
    n_c = formula.count("C")
    n_h = formula.count("H")
    if n_h != 2 * n_c:
        return False, f"Molecular formula {formula} does not match CnH2n"

    # Check for branching
    branched = False
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 3:
            branched = True
            break
    if not branched:
        return True, "Unbranched alkene"
    else:
        return True, "Branched alkene"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16236',
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.',
        'parents': ['CHEBI:33833', 'CHEBI:24631', 'CHEBI:51768']
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
    'num_true_positives': 350,
    'num_false_positives': 0,
    'num_true_negatives': 182114,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}