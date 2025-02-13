"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone skeleton SMARTS pattern
    flavanone_pattern = Chem.MolFromSmarts("c1ccccc1C2=CC(=O)CC(O2)c3ccccc3")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone skeleton found"
    
    # Define the 3'-hydroxy group SMARTS pattern on phenyl ring
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    matches = mol.GetSubstructMatches(hydroxy_3_prime_pattern)
    if not matches:
        return False, "No hydroxy group found at position 3' on the phenyl ring"
    
    return True, "Contains flavanone skeleton with a hydroxy group at the 3' position on the phenyl ring"