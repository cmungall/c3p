"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: CHEBI:28738 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid with a carboxyl group at C-24.
    It may contain glycine or taurine conjugates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 5beta-cholanic acid backbone
    cholanic_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@H]3[C@@H]4C[C@H](O)[C@@]5([C@@](CC[C@]5(C)C4)(C)C3)[C@@]12C)C")
    if not mol.HasSubstructMatch(cholanic_pattern):
        return False, "No 5beta-cholanic acid backbone found"
    
    # Check for carboxyl group at C-24
    carboxyl_pattern = Chem.MolFromSmarts("CCC(O)=O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, "Missing or multiple carboxyl groups at C-24"
    
    # Check for hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 2:
        return False, "Less than 2 hydroxy groups"
    
    # Check for glycine or taurine conjugates (optional)
    glycine_pattern = Chem.MolFromSmarts("NCC(O)=O")
    taurine_pattern = Chem.MolFromSmarts("S(=O)(=O)NCCC(O)=O")
    if mol.HasSubstructMatch(glycine_pattern) or mol.HasSubstructMatch(taurine_pattern):
        return True, "Contains 5beta-cholanic acid backbone with carboxyl group at C-24 and glycine/taurine conjugate"
    
    return True, "Contains 5beta-cholanic acid backbone with carboxyl group at C-24"