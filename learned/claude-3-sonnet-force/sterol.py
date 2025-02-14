"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:16829 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the tetracyclic steroid backbone
    backbone_pattern = Chem.MolFromSmarts("[C@@]12([C@]([C@]3([C@]([C@]1([C@@]4([C@](C[C@H](CC4)O)(CC3)[H])C)[H])[H])[H])(CC[C@@]([C@@H]([C@@H]2[H])C)[H])[H])C"
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No steroid backbone found"
    
    # Look for 3-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[O;H1]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, "Not exactly one hydroxy group found"
    
    # Check if hydroxy group is at the 3-position
    hydroxy_atom = mol.GetAtomWithIdx(hydroxy_matches[0][0])
    if hydroxy_atom.GetNeighbors()[0].GetAtomicNum() != 6 or hydroxy_atom.GetNeighbors()[1].GetAtomicNum() != 6:
        return False, "Hydroxy group not at 3-position"
    
    # Check for side chains
    side_chain_pattern = Chem.MolFromSmarts("[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) > 1:
        return True, "Contains steroid backbone with 3-hydroxy group and additional side chain(s)"
    
    return True, "Contains steroid backbone with 3-hydroxy group"