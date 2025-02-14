"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone of general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) ≠ H)
    in which the C=O function is conjugated to a C=C double bond at the alpha,beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the enone SMARTS pattern
    enone_pattern = Chem.MolFromSmarts('[C;H0;X3](=O)-[C]=[C]')
    if enone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha,beta-unsaturated ketone (enone) structure"
    else:
        return False, "No alpha,beta-unsaturated ketone (enone) structure found"