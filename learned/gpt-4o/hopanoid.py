"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is characterized by a hopane skeleton, which consists of five interconnected rings forming a specific triterpenoid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Hopane core: hopanoids are based on this C30 structure with 5 fused rings
    # Recognizing the basic polycyclic structure (triterpenoid characteristic),
    hopane_pattern = Chem.MolFromSmarts("C1[C@@H]2C[C@]34C(C2)C(C)(C3)CC4")  # This is a simplistic hopane pattern

    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains hopane skeleton"
    else:
        return False, "No hopane skeleton recognized"

# Example usage:
# print(is_hopanoid("CC(=C)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3C[C@H](O)[C@@]12C"))