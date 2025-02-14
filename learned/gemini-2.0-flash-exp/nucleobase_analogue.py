"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analog based on its SMILES string.
    A nucleobase analog is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analog, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the core pyrimidine or purine ring system
    pyrimidine_pattern = Chem.MolFromSmarts("n1cncc1")
    purine_pattern = Chem.MolFromSmarts("n1c2ncncn2c1")


    if not mol.HasSubstructMatch(pyrimidine_pattern) and not mol.HasSubstructMatch(purine_pattern):
         # 2. Check for specific Aza/Oxa analogs (check positions of N/O in the rings, next to a Carbon)
        azapurine_pattern = Chem.MolFromSmarts("n1[c,n]2[n][c][n][c]2[c,n]1")
        azapyrimidine_pattern = Chem.MolFromSmarts("n1[c][n][c][n][c]1")
        oxapurine_pattern = Chem.MolFromSmarts("o1[c,n]2[n][c][n][c]2[c,n]1")
        oxapyrimidine_pattern = Chem.MolFromSmarts("o1[c][n][c][n][c]1")
    
        if not mol.HasSubstructMatch(azapurine_pattern) and not mol.HasSubstructMatch(azapyrimidine_pattern) and not mol.HasSubstructMatch(oxapurine_pattern) and not mol.HasSubstructMatch(oxapyrimidine_pattern):
            return False, "No canonical pyrimidine or purine ring, or aza/oxa analog, found."

    return True, "Contains a pyrimidine or purine ring."