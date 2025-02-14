"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Lipopolysaccharides are complex molecules consisting of a trisaccharide repeating unit
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 
    3-hydroxytetradecanoic acid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool or None: True if molecule is a lipopolysaccharide, False otherwise or None if classification is not possible
        str or None: Reason for classification or None if not classified
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Example patterns; these would need refinement for genuine attempts
    # Looking for sugar alcohol structures and long carbon chains
    patterns = {
        "sugar": Chem.MolFromSmarts("C(C(C(C(C(C=O)O)O)O)O)O"),  # Simple sugar-like pattern
        "long_chain": Chem.MolFromSmarts("C(~C)~C~C~C~C~C~C~C~C"),  # Long chain carbon pattern
    }

    matched_patterns = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.append(name)

    if "sugar" in matched_patterns and "long_chain" in matched_patterns:
        # Simplistic condition to return True
        return True, "Contains indicative sugar and long-chain motifs commonly found in lipopolysaccharides"

    # Complexity likely surpassing this modality without specialized software/database access
    if not matched_patterns:
        return False, "No matching indicative patterns found. Complexity likely surpasses current SMILES-based method."

    return False, "Insufficient pattern match certainty for classification"