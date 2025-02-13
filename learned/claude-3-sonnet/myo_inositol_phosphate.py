"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons (should be exactly 6 for myo-inositol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Found {c_count} carbons, need exactly 6 for myo-inositol"
    
    # Check for phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count == 0:
        return False, "No phosphorus atoms found"
    
    # Check for ring size (should be 6)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size == 6 for size in ring_sizes):
        return False, "No 6-membered ring found"
        
    # Multiple patterns to match different representations of myo-inositol core
    myo_patterns = [
        # Pattern 1: Standard myo-inositol with explicit stereochemistry
        "[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
        # Pattern 2: Alternative representation
        "[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)1",
        # Pattern 3: More general pattern with any substituents
        "C1C(O)C(O)C(O)C(O)C(O)1"
    ]
    
    found_myo = False
    for pattern in myo_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_myo = True
            break
            
    if not found_myo:
        return False, "Does not match myo-inositol core pattern"
    
    # Check for phosphate groups - multiple patterns to match different forms
    phosphate_patterns = [
        # Standard phosphate
        "[OX2][P]([OX2,OX1-])([OX2,OX1-])=[OX1,OX2-]",
        # Charged phosphate
        "[O-][P](=[O])([O-])[O-]",
        # Alternative phosphate representation
        "[P](=O)([O,O-])([O,O-])[O,O-]"
    ]
    
    phosphate_found = False
    for pattern in phosphate_patterns:
        phos_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(phos_pattern):
            phosphate_found = True
            break
            
    if not phosphate_found:
        return False, "No phosphate groups found"
    
    # Verify phosphate attachment to ring
    ring_atoms = ring_info.AtomRings()[0]
    ring_carbons = set(idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    phosphate_attached = False
    for carbon_idx in ring_carbons:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetAtomicNum() == 15:  # Phosphorus
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
        if phosphate_attached:
            break
            
    if not phosphate_attached:
        return False, "Phosphate groups not properly attached to inositol ring"

    # Count total oxygens to verify complete phosphate groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    min_o_needed = 6 + (p_count * 3)  # 6 for ring hydroxyls + 3 per phosphate
    if o_count < min_o_needed:
        return False, "Insufficient oxygen atoms for complete phosphate groups"

    return True, f"Contains myo-inositol core with {p_count} phosphorus atoms in phosphate groups"