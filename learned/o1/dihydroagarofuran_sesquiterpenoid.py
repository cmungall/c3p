"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid (15-carbon terpenoid) with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the dihydroagarofuran skeleton as a SMARTS pattern
    # Note: This is a simplified pattern and may not capture all variations
    dihydroagarofuran_smarts = """
    [C@@H]1[C@H](C[C@@H]2[C@H]([C@H]1)O[C@H]3[C@@H](CC[C@@H]23)O)C
    """
    skeleton = Chem.MolFromSmarts(dihydroagarofuran_smarts)
    if skeleton is None:
        return False, "Invalid SMARTS pattern for dihydroagarofuran skeleton"

    # Check if the molecule contains the dihydroagarofuran skeleton
    if not mol.HasSubstructMatch(skeleton):
        return False, "Dihydroagarofuran skeleton not found"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Number of carbons ({num_carbons}) less than 15"

    # Check if it's a sesquiterpenoid (typically 15 carbons)
    if num_carbons >= 15:
        return True, "Contains dihydroagarofuran skeleton and has at least 15 carbons - classified as dihydroagarofuran sesquiterpenoid"
    else:
        return False, "Does not meet criteria for dihydroagarofuran sesquiterpenoid"