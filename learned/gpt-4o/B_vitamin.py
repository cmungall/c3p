"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Classifies a molecule as a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for B vitamins
    b_vitamins_patterns = {
        'Vitamin B1 (thiamine)': '[n+]-C-C-S',
        'Vitamin B2 (riboflavin)': 'C1=C(C=CC2=C1NC3=C(C=CC(=C3N=C2)O)O)O',
        'Vitamin B3 (niacin)': 'C1=CC=C(C=C1)C(=O)N',
        'Vitamin B5 (pantothenic acid)': 'C(C(C(C(=O)N)CO)C(=O)O)C',
        'Vitamin B6 (pyridoxine and derivatives)': 'c1c(CO)cnc(C)c1O',
        'Vitamin B7 (biotin)': '[H][C@@]1(C2=SCC[NH2+]2)[C@@H](C(=O)NC1)CCCCC(=O)O',
        'Vitamin B9 (folic acid and derivatives)': 'NC1=C(NC=NC2=C1N=C(NC2=O)N)C3=CC=CC=C3',
        'Vitamin B12 (cobalamin-related structures)': 'CNC1=CC2=C(C=C1)N=C(N=O2)C=3[C@H]([C@H]([C@H](O3)O)O)C'
    }
    
    # Check for substructures
    for name, pattern in b_vitamins_patterns.items():
        substruct = Chem.MolFromSmarts(pattern)
        if not substruct:
            continue
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches substructure for {name}"

    return False, "No B vitamin substructure detected"