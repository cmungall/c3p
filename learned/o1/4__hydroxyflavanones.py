"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a flavanone with a hydroxy substituent at the 4' position of the B-ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS pattern
    flavanone_core_smarts = 'C1(=O)CCOc2ccccc12'  # Flavanone core with rings fused at positions 1 and 2
    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if flavanone_core is None:
        return None, "Error in flavanone core SMARTS pattern"

    # Check if the molecule contains the flavanone core
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Flavanone core not found"

    # Define SMARTS pattern for hydroxy group at 4' position on B-ring
    hydroxy_b_ring_smarts = 'C1(=O)CCOc2ccc(O)cc12'  # Hydroxy at para position on B-ring
    hydroxy_b_ring = Chem.MolFromSmarts(hydroxy_b_ring_smarts)
    if hydroxy_b_ring is None:
        return None, "Error in 4'-hydroxyflavanone SMARTS pattern"

    # Check if molecule matches the 4'-hydroxyflavanone pattern
    if mol.HasSubstructMatch(hydroxy_b_ring):
        return True, "Molecule is a 4'-hydroxyflavanone"
    else:
        return False, "4'-hydroxy group at B-ring not found"