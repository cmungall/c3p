"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic flavonoid core patterns
    flavone_pattern = Chem.MolFromSmarts("[#6]1=[#6][#6](=[O])[#6]2=[#6][#6]=[#6][#6]=[#6]2[O][#6]1")  # C-C=C(=O)-C1=CC=CC=C1-O-C
    flavan_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2=[#6][#6]=[#6][#6]=[#6]2[O][#6]1")  # Basic flavan core
    
    # Count flavonoid cores
    flavone_matches = len(mol.GetSubstructMatches(flavone_pattern))
    flavan_matches = len(mol.GetSubstructMatches(flavan_pattern))
    total_cores = flavone_matches + flavan_matches
    
    if total_cores < 2:
        return False, f"Found only {total_cores} flavonoid cores, need at least 2"
    
    # Check for characteristic oxygen atoms (typically from hydroxyl/methoxy groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:  # Minimum oxygens expected in a biflavonoid
        return False, f"Too few oxygen atoms ({o_count}) for a biflavonoid"
    
    # Check for aromatic rings (should have multiple)
    aromatic_rings = 0
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
    
    if aromatic_rings < 4:  # Biflavonoids typically have at least 4 aromatic rings
        return False, f"Too few aromatic rings ({aromatic_rings}) for a biflavonoid"
    
    # Check molecular weight (biflavonoids are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:  # Typical biflavonoids are >450 Da
        return False, "Molecular weight too low for biflavonoid"
        
    # Check for hydroxyl groups (most biflavonoids have multiple OH groups)
    oh_pattern = Chem.MolFromSmarts("[OH]")
    oh_count = len(mol.GetSubstructMatches(oh_pattern))
    if oh_count < 2:
        return False, f"Too few hydroxyl groups ({oh_count}) for a biflavonoid"
    
    # Check carbon count (biflavonoids typically have 30-40 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 45:
        return False, f"Carbon count ({c_count}) outside typical range for biflavonoids"
    
    # If all checks pass, it's likely a biflavonoid
    return True, f"Contains {total_cores} flavonoid cores with appropriate substitution pattern"