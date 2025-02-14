"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:35620 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyanohydrin substructure: R-C(OH)(C#N)-R'
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4]([OH])([CX2]#N)")
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "Missing cyanohydrin substructure"
    
    # Check if the cyanohydrin substructure is formed by addition to a carbonyl compound
    cyanohydrin_atoms = mol.GetSubstructMatches(cyanohydrin_pattern)[0]
    carbon_atom = mol.GetAtomWithIdx(cyanohydrin_atoms[0])
    
    # Look for patterns indicating addition to a carbonyl
    carbonyl_pattern1 = Chem.MolFromSmarts("[CX3]([CX4]([OH])([CX2]#N))")  # R-C(=O)-C(OH)(C#N)-R'
    carbonyl_pattern2 = Chem.MolFromSmarts("[CX3H1]([CX4]([OH])([CX2]#N))")  # R-CH(=O)-C(OH)(C#N)-R'
    
    if mol.HasSubstructMatch(carbonyl_pattern1) or mol.HasSubstructMatch(carbonyl_pattern2):
        return True, "Contains the cyanohydrin substructure formed by addition of hydrogen cyanide to a carbonyl group"
    
    # If no pattern indicating addition to a carbonyl is found, assume it's a cyanohydrin
    return True, "Contains the cyanohydrin substructure"