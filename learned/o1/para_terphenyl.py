"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl consists of three benzene rings connected linearly at the
    para positions (positions 1 and 4) of the central ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define para-terphenyl core SMARTS pattern
    para_terphenyl_smarts = "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"
    pattern = Chem.MolFromSmarts(para_terphenyl_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Contains para-terphenyl core structure"
    else:
        return False, "Does not contain para-terphenyl core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'para-terphenyl',
        'definition': 'A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.'
    }
}