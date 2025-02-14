"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:24351 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import Counter

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion has an oxo group at the 2-position and a negatively charged carboxylate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group [-C(=O)[O-]]
    carboxylate_pattern = Chem.MolFromSmarts("[C-](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, need exactly 1"
    
    # Look for oxo group [C(=O)-]
    oxo_pattern = Chem.MolFromSmarts("[C](=O)(-*)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) != 1:
        return False, f"Found {len(oxo_matches)} oxo groups, need exactly 1"
    
    # Check if oxo group is at 2-position
    carboxylate_atom = mol.GetAtomWithIdx(carboxylate_matches[0][0])
    oxo_atom = mol.GetAtomWithIdx(oxo_matches[0][0])
    
    if carboxylate_atom.GetNeighbors()[0].GetIdx() != oxo_atom.GetIdx():
        return False, "Oxo group not at 2-position relative to carboxylate"
    
    # Check for only one carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[C](=O)[O;!-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) > 0:
        return False, "Contains non-ionized carboxylic acid groups"
    
    # Count atoms to ensure monocarboxylic
    atom_counts = Counter(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if atom_counts[6] > 10:  # Allow up to 10 carbon atoms
        return False, "Too many carbon atoms for monocarboxylic acid"
    
    return True, "Contains a carboxylate group and an oxo group at the 2-position"