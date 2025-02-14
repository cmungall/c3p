"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains an epoxide ring and a long aliphatic chain ending in a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an epoxide group (cyclic ether C1OC1)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")  # General cyclic ether for epoxides
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"
    
    # Check for the carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a sufficient aliphatic chain length (at least 10 carbons)
    aliphatic_carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() <= 4)
    if aliphatic_carbon_count < 10:
        return False, "Aliphatic chain too short to be a fatty acid"

    return True, "Contains epoxide ring as part of a long aliphatic chain ending in a carboxyl group"