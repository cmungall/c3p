"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: CHEBI:27913 steroid sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the formal condensation
    of a hydroxy group of any steroid with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@]1(CC[C@@]2([H])[C@]3([H])C[C@H]([C@H]4[C@@]5([C@H]([C@@H]([C@@H]5CC=C6C[C@@H]4C=C7[C@@H]6C[C@@H]7[H])O)[H])C)=C3CCC2=O)C[C@@H]1[H]"
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
    
    # Look for sulfate group pattern (-OS(=O)(=O)O)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[OX2]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate group found"
    
    # Check if sulfate group is attached to an oxygen of a hydroxy group
    for match in sulfate_matches:
        sulfate_atom = mol.GetAtomWithIdx(match[0])
        sulfate_neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in sulfate_atom.GetNeighbors()]
        if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1 for neighbor in sulfate_neighbors):
            # Sulfate group is attached to a hydroxy group
            return True, "Contains a steroid backbone with a sulfate group attached to a hydroxy group"
    
    return False, "Sulfate group is not attached to a hydroxy group"