"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Key molecular patterns
    sphingoid_base_pattern = Chem.MolFromSmarts("N[C@@H](CO)C([C@H](O)C(O)C)")  # General pattern for long-chain amino alcohol
    amide_pattern = Chem.MolFromSmarts("C(=O)N")  # General amide linkage pattern
    
    # Match patterns
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base pattern found"
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Check for sufficient carbon content (indicative of the long-chain nature of ceramides)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:
        return False, f"Too few carbon atoms ({c_count}), typical for ceramides"
    
    # Flexible functionality to identify variations like hydroxyl groups
    hydroxyl_group_pattern = Chem.MolFromSmarts("[OX2H]")  # Smart pattern to match -OH groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_group_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found, common in ceramides"
    
    return True, "Contains characteristic sphingoid base with an amide-linked fatty acid chain"