"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:33716 dihydroflavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroflavonol(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the flavanone scaffold
    flavanone_pattern = Chem.MolFromSmarts("[O;R]1C(=O)C2=C(C=CC=C2)C2=C1C=CC(O)=C2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Molecule does not contain the flavanone scaffold"
    
    # Check for hydroxy group at position 3 of the heterocyclic ring
    dihydroflavonol_pattern = Chem.MolFromSmarts("[O;R]1C(=O)C2=C(C=CC=C2)C2=C1C(O)=CC=C2")
    if not mol.HasSubstructMatch(dihydroflavonol_pattern):
        return False, "No hydroxy group at position 3 of the heterocyclic ring"
    
    return True, "Molecule contains the flavanone scaffold with a hydroxy group at position 3 of the heterocyclic ring"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33716',
        'name': 'dihydroflavonol',
        'definition': 'Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.',
        'parents': ['CHEBI:28519', 'CHEBI:35402']
    },
    # Additional metadata can be added here
}