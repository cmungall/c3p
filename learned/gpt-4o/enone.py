"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the general formula:
    R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) != H), where the C=O is conjugated to a C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Improved SMARTS for enone: broadening conditions and patterns.
    # Each atom is specified with constraints: 
    # - Two carbon double bond (C=C)
    # - Followed by a single bond ketone (C=O) with non-H
    enone_pattern = Chem.MolFromSmarts("[C;R0]=[C;R0]-[C;R0](=O)[C,N,#6]")

    # Check for enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha, beta-unsaturated ketone (enone) structure"
    
    return False, "Does not contain alpha, beta-unsaturated ketone (enone) structure"