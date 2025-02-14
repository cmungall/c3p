"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the general formula R1R2C=CR3-C(=O)R4 (R4 != H)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the enone motif
    # The pattern now specifies that the carbonyl carbon is not bonded to any H (except via double bond)
    enone_pattern = Chem.MolFromSmarts("[C;!H0]=[C;!H0]-[C](=[O])[C;!H0]")
    
    # Find the enone substructure
    matches = mol.GetSubstructMatches(enone_pattern)
    
    if matches:
        return True, "Contains an alpha,beta-unsaturated ketone (C=C-C=O) motif with R4 != H"
    
    return False, "No alpha,beta-unsaturated ketone (C=C-C=O) motif found."