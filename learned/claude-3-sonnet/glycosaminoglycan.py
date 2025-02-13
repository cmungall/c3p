"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:36975 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is any polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aminomonosaccharide residues
    aminosugar_pattern = Chem.MolFromSmarts("[NX3H2][CR]")
    aminosugar_matches = mol.GetSubstructMatches(aminosugar_pattern)
    
    # Check for polysaccharide backbone
    saccharide_pattern = Chem.MolFromSmarts("[OX2]C[OX2]")
    saccharide_matches = mol.GetSubstructMatches(saccharide_pattern)
    
    # Classify as glycosaminoglycan if both patterns are found
    if aminosugar_matches and saccharide_matches:
        return True, "Contains aminomonosaccharide residues in a polysaccharide backbone"
    else:
        return False, "Does not contain aminomonosaccharide residues in a polysaccharide backbone"