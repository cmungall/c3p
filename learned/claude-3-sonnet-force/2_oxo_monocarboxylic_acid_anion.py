"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:24351 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion has an oxo group and a negatively charged carboxylate group,
    with the oxo group being in the 2-position relative to the carboxylate group.

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
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, need exactly 1"
    
    # Look for oxo group [C(=O)-]
    oxo_pattern = Chem.MolFromSmarts("[C](=O)(-*)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "No oxo groups found"
    
    carboxylate_atom = mol.GetAtomWithIdx(carboxylate_matches[0][0])
    
    # Check if any oxo group is in the 2-position relative to the carboxylate
    for oxo_match in oxo_matches:
        oxo_atom = mol.GetAtomWithIdx(oxo_match[0])
        path = Chem.rdmolops.GetShortestPath(mol, carboxylate_atom.GetIdx(), oxo_atom.GetIdx())
        if len(path) >= 3:
            break
    else:
        return False, "No oxo group found in the 2-position relative to the carboxylate"
    
    # Check for multiple carboxylate groups
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) > 1:
        return False, "Multiple carboxylate groups found, should be monocarboxylic"
    
    return True, "Contains a carboxylate group and an oxo group in the 2-position"