"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a phenyl-substituted benzopyran structure, often with C15 or C16 skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is likely a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for various flavonoid subclasses
    flavonoid_core_patterns = [
        # Basic flavonoid structure
        Chem.MolFromSmarts("c1cc(O)c2c(c1)C=CC(=O)O2"),
        # Isoflavonoid structure
        Chem.MolFromSmarts("c1cc2c(cc1)O[C]=C(C=O)C2"),
        # Chalcone structure
        Chem.MolFromSmarts("c1cc(O)cc(c1)C=CC(=O)c2ccccc2"),
        # Neoflavonoid structure (3-arylchromen-4-one)
        Chem.MolFromSmarts("c1ccc2c(c1)C=CC(=O)O2"),
        # Flavanone
        Chem.MolFromSmarts("c1ccc2c(c1)OC(=O)C(C2)O"),
    ]
    
    # Search for flavonoid core structures
    for pattern in flavonoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid substructure (phenyl-substituted chromanone or related pattern)"

    return False, "No core flavonoid substructure found"