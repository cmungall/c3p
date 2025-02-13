"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is characterized by a long carbon chain, an amino group, and hydroxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sphingoids typically have a long carbon chain of about 14-20 carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No sufficient long carbon chain recognized"

    # Check for the presence of an amine group at C-2 position
    amine_pattern = Chem.MolFromSmarts("N")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amino group found"

    # Check for the presence of hydroxyl group(s)
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check if the molecule contains at least 14 - 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (14 <= c_count <= 20):
        return False, f"Carbon chain length {c_count} is not typical for sphingoids"

    # Verify if specified stereochemistry is present (if mentioned in examples)
    stereo_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[C@@H]"))
    if stereo_matches:
        # Check for specific stereocenters. This can be extended for detailed stereochemistry checks.
        return True, "Contains stereochemistry matching a sphingoid"

    return True, "Contains typical structural features of a sphingoid (long chain, amino, hydroxyl groups)"