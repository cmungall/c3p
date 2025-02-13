"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide contains a macrocyclic lactone with a ring of twelve or more members.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify large rings (12 or more members)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    has_large_ring = any(size >= 12 for size in ring_sizes)
    if not has_large_ring:
        return False, "No ring with 12 or more members found"

    # Look for lactone (cyclic ester) presence
    lactone_smarts = "[OX2H1][#6]1~[#6](=[OX1])[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1" # 12-member ring with an ester
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No macrocyclic lactone structure found"
    
    return True, "Contains a macrocyclic lactone ring with 12 or more members"