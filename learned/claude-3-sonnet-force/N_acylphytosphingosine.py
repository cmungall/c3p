"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:35801 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phytosphingosine backbone pattern
    phytosphingosine_pattern = Chem.MolFromSmarts("[C@H](CO)[C@@H](O)[C@@H](O)CCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Look for acyl group attached to nitrogen
    acyl_pattern = Chem.MolFromSmarts("[N;H1][C;X3](=[O;X1])[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to nitrogen found"
    
    # Look for galactose or glucose moiety
    sugar_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No galactose or glucose moiety found"
    
    return True, "Contains phytosphingosine backbone with acyl group attached to nitrogen and galactose/glucose moiety"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35801',
        'name': 'N-acylphytosphingosine',
        'definition': 'A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.',
        'parents': ['CHEBI:35798']
    },
    'config': {
        # ... (existing config data)
    },
    # ... (existing metadata)
}