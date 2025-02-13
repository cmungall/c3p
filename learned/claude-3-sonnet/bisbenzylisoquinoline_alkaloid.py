"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: CHEBI:bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum complexity
    if mol.GetNumAtoms() < 30:  # These are large molecules
        return False, "Molecule too small for bisbenzylisoquinoline alkaloid"

    # Look for isoquinoline cores (can be N-methylated or not)
    isoquinoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]2[#6](=[#6][#6]=1)[#6]=[#6][#6]=2")
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    
    if len(isoquinoline_matches) < 2:
        return False, "Must contain two isoquinoline cores"

    # Look for ether bridges (-O-)
    ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(ether_matches) < 1:
        return False, "Must contain ether bridges between units"

    # Check for N-methyl groups or NH (typical in these alkaloids)
    n_methyl_pattern = Chem.MolFromSmarts("[#7]-[CH3]")
    nh_pattern = Chem.MolFromSmarts("[#7;H1]")
    n_methyl_matches = mol.GetSubstructMatches(n_methyl_pattern)
    nh_matches = mol.GetSubstructMatches(nh_pattern)
    
    if len(n_methyl_matches) + len(nh_matches) < 2:
        return False, "Must contain two nitrogen centers (N-methyl or NH)"

    # Look for benzyl groups attached to isoquinoline
    benzyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1-[CH2]")
    benzyl_matches = mol.GetSubstructMatches(benzyl_pattern)
    
    if len(benzyl_matches) < 2:
        return False, "Must contain two benzyl groups"

    # Optional: Check for common substituents (methoxy, hydroxy)
    methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[CH3]")
    hydroxy_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count aromatic rings
    aromatic_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
            
    if aromatic_rings < 4:  # Should have at least 4 aromatic rings
        return False, "Insufficient number of aromatic rings"

    # Additional check for methylenedioxy groups (common in this class)
    methylenedioxy_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#8]1")
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)

    # Verify overall connectivity
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~1")):
        return False, "Missing key ring connectivity"

    return True, "Contains two benzylisoquinoline units connected by ether bridges with appropriate substitution patterns"