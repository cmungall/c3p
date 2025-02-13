"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    This involves checking for more than ten glycosidic bonds (O-C-O structure in sugar rings) 
    and a polysaccharide backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar ring pattern (typically pyranose or furanose forms)
    # A simple ring count for demonstration purposes
    sugar_ring_pattern = Chem.MolFromSmarts("O1C(O)C(O)C(O)C(O)C1")
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)

    # We assume these matches represent a basic sugar unit,
    # though polysaccharides might have complex patterns
    if len(sugar_matches) <= 10:
        return False, f"Found {len(sugar_matches)} monosaccharide rings, need more than 10"

    # Check for connectivity between rings - needs adjustment based on real-world polysaccharide features
    # For simplification: Check for overall number of ethers (O-C-O types)
    ether_pattern = Chem.MolFromSmarts("[OX2]C[OX2]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(ether_matches) <= 10:
        return False, f"Found {len(ether_matches)} glycosidic bonds, need more than 10"

    return True, "Contains more than ten monosaccharide rings linked via glycosidic bonds"