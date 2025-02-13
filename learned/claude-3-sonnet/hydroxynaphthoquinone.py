"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthaoquinone moiety is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_hydroxynaphthoquinone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for naphthoquinone core
    # Pattern matches two fused 6-membered rings with two ketone groups
    naphthoquinone_pattern = Chem.MolFromSmarts("[#6]1[#6]=[#6][#6](=[O])[#6]2[#6][#6][#6][#6][#6]2[#6]1=[O]")
    
    # Alternative pattern for different ketone arrangements
    naphthoquinone_pattern2 = Chem.MolFromSmarts("[#6]1[#6](=[O])[#6]=[#6][#6]2[#6][#6][#6][#6][#6]2[#6]1=[O]")
    
    if not (mol.HasSubstructMatch(naphthoquinone_pattern) or mol.HasSubstructMatch(naphthoquinone_pattern2)):
        return False, "No naphthoquinone core found"
    
    # Check for hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if not hydroxy_matches:
        return False, "No hydroxy groups found"
    
    # Get the naphthoquinone core atoms
    core_matches = mol.GetSubstructMatches(naphthoquinone_pattern)
    if not core_matches:
        core_matches = mol.GetSubstructMatches(naphthoquinone_pattern2)
    core_atoms = set(core_matches[0])
    
    # Check if any hydroxy group is attached to the naphthoquinone core
    for hydroxy_match in hydroxy_matches:
        oh_oxygen = hydroxy_match[0]
        # Get the atom the OH is attached to
        neighbors = mol.GetAtomWithIdx(oh_oxygen).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() in core_atoms:
                return True, "Contains naphthoquinone core with at least one hydroxy group attached"
    
    return False, "Hydroxy group(s) present but not attached to naphthoquinone core"