"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:36314 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide is a natural compound consisting of:
    - A trisaccharide repeating unit (two heptose units and octulosonic acid)
    - Oligosaccharide side chains
    - 3-hydroxytetradecanoic acid units
    They are a major constituent of the cell walls of Gram-negative bacteria.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for trisaccharide repeating unit
    heptose_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@H]([C@@H]([C@@H](C1)O)O)O)O")
    octulosonic_acid_pattern = Chem.MolFromSmarts("O=C(O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO")
    trisaccharide_matches = mol.GetSubstructMatches(heptose_pattern) + mol.GetSubstructMatches(octulosonic_acid_pattern)
    if len(trisaccharide_matches) < 3:
        return False, f"Missing trisaccharide repeating unit, got {len(trisaccharide_matches)} matches"
    
    # Look for oligosaccharide side chains
    oligosaccharide_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    oligosaccharide_matches = mol.GetSubstructMatches(oligosaccharide_pattern)
    if len(oligosaccharide_matches) < 3:
        return False, f"Missing oligosaccharide side chains, got {len(oligosaccharide_matches)} matches"
    
    # Look for 3-hydroxytetradecanoic acid unit
    hydroxy_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC[C@@H](O)C(O)=O")
    hydroxy_acid_matches = mol.GetSubstructMatches(hydroxy_acid_pattern)
    if len(hydroxy_acid_matches) == 0:
        return False, "Missing 3-hydroxytetradecanoic acid unit"
    
    # Check molecular weight - LPS typically >2000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 2000:
        return False, "Molecular weight too low for lipopolysaccharide"
    
    return True, "Contains trisaccharide repeating unit, oligosaccharide side chains, and 3-hydroxytetradecanoic acid units"