"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: CHEBI:35468 anthoxanthin

Anthoxanthins are a type of flavonoid pigments in plants. They are water-soluble pigments
which range in color from white or colorless to a creamy to yellow, often on petals of flowers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is a flavonoid
    flavonoid_pattern = Chem.MolFromSmarts("c1c(c2c(cc1)occc2)ccc(O)c(O)c")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Not a flavonoid structure"
    
    # Check for common anthoxanthin functional groups/substituents
    hydroxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    methoxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRingSize(2))
    sulfate_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16 and atom.GetTotalNumHs() == 0)
    glucoside_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRingSize(6) and atom.GetTotalNumHs() == 1)
    
    # Common anthoxanthin substituents: hydroxy, methoxy, sulfate, glucoside
    if hydroxy_count >= 3 or methoxy_count >= 1 or sulfate_count >= 1 or glucoside_count >= 1:
        return True, "Contains common anthoxanthin functional groups/substituents"
    
    return False, "Does not contain typical anthoxanthin features"