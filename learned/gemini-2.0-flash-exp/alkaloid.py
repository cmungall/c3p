"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is likely an alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen in a ring
    nitrogen_in_ring_pattern = Chem.MolFromSmarts("[#7;R]")
    if not mol.HasSubstructMatch(nitrogen_in_ring_pattern):
        return False, "No nitrogen in a ring"

    # Exocyclic nitrogen check. This is a bit tricky as we do not want amines
    exocyclic_nitrogen_pattern = Chem.MolFromSmarts("[#7;!R]")
    if mol.HasSubstructMatch(exocyclic_nitrogen_pattern):
        return False, "Contains exocyclic nitrogen"
    
    # Check for carboxylic acid group directly linked to nitrogen in a ring (amino acid check)
    amino_acid_pattern = Chem.MolFromSmarts("[#7;R]-[C](=[O])[OH]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Contains carboxylic acid group directly linked to nitrogen in a ring"

    # Check for carboxylic acid group linked to a carbon with a nitrogen attached
    amino_acid_pattern2 = Chem.MolFromSmarts("[#7;R]-[C]-[C](=[O])[OH]")
    if mol.HasSubstructMatch(amino_acid_pattern2):
        return False, "Contains carboxylic acid group linked to a carbon with a nitrogen attached"
        
    return True, "Contains at least one nitrogen in a ring, does not have exocyclic nitrogen or carboxylic acid groups"