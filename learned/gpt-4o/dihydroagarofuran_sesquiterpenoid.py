"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is any sesquiterpenoid with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a dihydroagarofuran skeleton
    # This SMARTS is a simplified representation; real implementations might require fine-tuning
    dihydroagarofuran_pattern = Chem.MolFromSmarts("C1(C)OC2CCC3CCC(C3)C(C2)O1")

    if not mol.HasSubstructMatch(dihydroagarofuran_pattern):
        return False, "No dihydroagarofuran skeleton found in molecule"

    # Further checks specific to the class may be added here if necessary

    return True, "Contains dihydroagarofuran skeleton"