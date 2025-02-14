"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: CHEBI:51841 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid containing a dihydroagarofuran skeleton.

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

    # Check if molecule is a sesquiterpenoid (C15H24 or C15H22)
    formula = AllChem.CalcMolFormula(mol)
    if formula not in ['C15H24', 'C15H22']:
        return False, "Not a sesquiterpenoid (formula not C15H24 or C15H22)"

    # Look for dihydroagarofuran skeleton pattern
    dihydroagarofuran_pattern = Chem.MolFromSmarts("[C@H]1[C@]2([C@@H](OC(=O)[C])OC(=O)[C])C[C@@H]([C@H](OC(=O)[C])[C@@H]1OC(=O)[C])[C@@H]2OC(=O)[C]")
    if not mol.HasSubstructMatch(dihydroagarofuran_pattern):
        return False, "No dihydroagarofuran skeleton found"

    return True, "Contains the dihydroagarofuran skeleton characteristic of dihydroagarofuran sesquiterpenoids"