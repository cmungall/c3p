"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a C28 steroid with a modified side chain forming a lactone ring, and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid scaffold (cyclopentanoperhydrophenanthrene, C27H45)
    # Using a simplified SMARTS pattern for identifying the steroid backbone
    steroid_pattern = Chem.MolFromSmarts('[#6]12CC[C@H]3[C@@H](C2)CCC4=C1C=CC(C4)=O')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for lactone ring presence, e.g., a basic 5-membered ring lactone
    lactone_pattern = Chem.MolFromSmarts('O=C1OC2CC1[CX3](=O)')
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Ensure the molecule is a C28 steroid (approximately)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 28:
        return False, f"Too few carbons for a C28 steroid, found {carbon_count}"
    
    return True, "Molecule has a steroid backbone, a lactone ring, and a C28 framework"