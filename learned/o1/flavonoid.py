"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is based on 1-benzopyran with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the flavonoid core SMARTS pattern
    # This pattern represents a benzopyran fused ring system with a phenyl ring attached at position 2
    flavonoid_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1            # A-ring: benzene ring
    -[#6]2:[#8]:[#6]:[#6]:[#6]:[#6]:2           # C-ring: pyran ring fused to A-ring
    -[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3           # B-ring: phenyl ring attached at position 2
    """
    # Remove whitespace and newlines from the SMARTS string
    flavonoid_smarts = "".join(flavonoid_smarts.split())
    
    # Create the SMARTS pattern
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)
    if flavonoid_pattern is None:
        return False, "Failed to construct flavonoid SMARTS pattern"
    
    # Check if the molecule matches the flavonoid core
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Molecule contains flavonoid core structure"
    else:
        return False, "Molecule does not contain flavonoid core structure"