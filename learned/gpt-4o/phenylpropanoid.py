"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Includes molecules derived from phenylpropane with aromatic rings and characteristic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect aromatic systems which are central to phenylpropanoids
    aromatic_system_pattern = Chem.MolFromSmarts("a1aaaaa1")
    if not mol.HasSubstructMatch(aromatic_system_pattern):
        return False, "No aromatic system found"

    # Check for phenylpropanoid-like linkages or C6-C3 structure
    phenylpropanoid_pattern = Chem.MolFromSmarts("c1ccccc1CC=C")  # Simplified phenylpropane topology
    if not mol.HasSubstructMatch(phenylpropanoid_pattern):
        return False, "No phenylpropanoid structure found"

    # Check for frequently occurring functional groups in phenylpropanoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # -OH groups
    ether_pattern = Chem.MolFromSmarts("CO")  # Methoxy groups
    carbonyl_pattern = Chem.MolFromSmarts("C=O")  # Ketones, aldehydes
    olefinic_pattern = Chem.MolFromSmarts("C=C")  # Double-bond in side chains

    has_functional_group = any(
        mol.HasSubstructMatch(pattern) 
        for pattern in [hydroxyl_pattern, ether_pattern, carbonyl_pattern, olefinic_pattern]
    )
    
    if not has_functional_group:
        return False, "Missing characteristic functional groups for phenylpropanoids"

    return True, "Contains aromatic system and expected phenylpropanoid linkages & functional groups"