"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by the intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group, resulting in a 1-oxacycloalkan-2-ol or unsaturated analogue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for lactol
    # Lactol functional group: a ring carbon bonded to an OH group and an ether oxygen in the ring
    lactol_smarts = '[C;R]([O;H1])[O;R]'
    lactol_pattern = Chem.MolFromSmarts(lactol_smarts)
    if lactol_pattern is None:
        return False, "Invalid lactol SMARTS pattern"

    # Search for lactol pattern in the molecule
    matches = mol.GetSubstructMatches(lactol_pattern)
    if matches:
        return True, "Contains lactol functional group"
    else:
        return False, "Lactol functional group not found"