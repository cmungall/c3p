"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid skeleton pattern (cyclopenta[a]phenanthrene) using chiral centers as proxies
    steroid_skeleton_pattern = Chem.MolFromSmarts("[C@]1([C@H](C2[C@H](C[C@H]3[C@@H](C[C@@H]4[C@]3(CO)C)C)C[C@]2(C)C)C)[C@H](O)C(C)CCC(=O)N")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroid skeleton typical of bile acids found"

    # Conjugation moiety patterns
    # Glycine pattern: amino acid with terminal COOH
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    # Taurine pattern: sulfonic acid
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)O")
    # Additional patterns can be added for other conjugate molecules (glucuronic acid, sugars, etc.)

    # Check for the presence of conjugation patterns
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Steroid skeleton with glycine conjugate found"
    if mol.HasSubstructMatch(taurine_pattern):
        return True, "Steroid skeleton with taurine conjugate found"

    # Check for sulfuric acid conjugates or sugar conjugates
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)")
    sugar_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H]([O-])C[C@@H]([O-])C([O-])SC1=O")
    if mol.HasSubstructMatch(sulfate_pattern):
        return True, "Steroid skeleton with sulfate conjugate found"
    if mol.HasSubstructMatch(sugar_pattern):
        return True, "Steroid skeleton with sugar conjugate found"

    return False, "No recognizable bile acid conjugate structure found"

__metadata__ = {   'chemical_class': {   'id': 'custom',
                          'name': 'bile acid conjugate',
                          'definition': 'Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.'},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}