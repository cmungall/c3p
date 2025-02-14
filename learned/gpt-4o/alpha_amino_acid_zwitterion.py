"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Defining a core zwitterionic pattern: amino group and carboxylate
    core_zwitterion_pattern = Chem.MolFromSmarts('[C@@H]([NH3+])C(=O)[O-]')
    
    # Expanded zwitterion patterns considering various complex amino acid structures
    extended_patterns = [
        Chem.MolFromSmarts('[NH3+]-[C@H](-C(=O)[O-])'),
        Chem.MolFromSmarts('[C@H]([NH3+])(C(=O)[O-])'),   # Alpha carbon with zwitterion
        Chem.MolFromSmarts('(C1O[C@@H](CN)[C@H](O1))C(=O)[O-]'),  # Ring structures (e.g., proline)
        Chem.MolFromSmarts('(N-=[N+]=N)-C')   # Attached beyond basic amine (extended chains)
    ]

    # First, check if the core structure is present
    if not mol.HasSubstructMatch(core_zwitterion_pattern):
        return False, "No matching alpha-amino acid zwitterion core structure found"
    
    # Check each stored pattern for additional confirmation of species complexity
    # Mark as zwitterion if any matching pattern is present
    for pattern in extended_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains alpha-amino acid structure with zwitterionic charges"
    
    return True, "Detected basic zwitterionic structure"

# Example usage:
# smiles = "CC(=O)CC([NH3+])C([O-])=O"  # Example SMILES from the list
# result, reason = is_alpha_amino_acid_zwitterion(smiles)
# print(f"Result: {result}, Reason: {reason}")