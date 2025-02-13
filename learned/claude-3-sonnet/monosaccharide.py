"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:16646 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with 3 or more carbon atoms
    and a (potential) carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count atoms and check for aldehyde/ketone group
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    has_carbonyl = any(atom.GetHybridization() == Chem.HybridizationType.SP2 and
                       sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2
                       for atom in mol.GetAtoms())
    
    # Check for 3 or more carbons, multiple oxygens, and carbonyl group
    if c_count < 3:
        return False, "Fewer than 3 carbon atoms"
    if o_count < 2:
        return False, "Fewer than 2 oxygen atoms"
    if not has_carbonyl:
        return False, "No carbonyl group found"
    
    # Look for polyhydroxy pattern (multiple -CH(O)-) and single carbonyl group
    polyhydroxy_pattern = Chem.MolFromSmarts("[CH1](O)")
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    polyhydroxy_matches = mol.GetSubstructMatches(polyhydroxy_pattern)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    if len(polyhydroxy_matches) < 2:
        return False, "Not enough hydroxyl groups"
    if len(carbonyl_matches) != 1:
        return False, "Incorrect number of carbonyl groups"
    
    # Check for cyclic monosaccharides
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Check for pyranose or furanose ring
        pyranose_pattern = Chem.MolFromSmarts("OC1CCCCO1")
        furanose_pattern = Chem.MolFromSmarts("OC1CCOC1")
        if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
            return True, "Cyclic monosaccharide detected"
        else:
            return False, "Unrecognized ring system"
    
    return True, "Acyclic monosaccharide detected"