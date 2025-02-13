"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:38047 1,2,4-triazine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine ring structure.
    1,2,4-triazine has nitrogen atoms at positions 1, 2, and 4 of a 6-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 1,2,4-triazine ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for 1,2,4-triazine core with explicit positions
    # Pattern 1: Aromatic form
    pattern1 = Chem.MolFromSmarts("[n]1[n]c[n]cc1")
    # Pattern 2: Non-aromatic form with double bonds
    pattern2 = Chem.MolFromSmarts("[N]1=[N][C]=[N][C]=[C]1")
    # Pattern 3: Alternative non-aromatic form
    pattern3 = Chem.MolFromSmarts("[N]1[N]=[C][N]=[C][C]1")
    
    patterns = [pattern1, pattern2, pattern3]
    match_found = False
    
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # Get the matched atoms
                ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
                
                # Check if this is part of a larger fused system
                is_fused = False
                for atom in ring_atoms:
                    # Count number of rings this atom belongs to
                    if len(mol.GetRingInfo().AtomRings()) > 1:
                        rings_count = 0
                        for ring in mol.GetRingInfo().AtomRings():
                            if atom.GetIdx() in ring:
                                rings_count += 1
                        if rings_count > 1:
                            is_fused = True
                            break
                
                # Verify nitrogen positions (1,2,4)
                n_positions = []
                for i, atom in enumerate(ring_atoms):
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        n_positions.append(i)
                
                # Check if nitrogens are in positions 0,1,3 (equivalent to 1,2,4)
                if sorted(n_positions) == [0,1,3] and not is_fused:
                    match_found = True
                    is_aromatic = all(atom.GetIsAromatic() for atom in ring_atoms)
                    ring_type = "aromatic" if is_aromatic else "non-aromatic"
                    return True, f"Contains {ring_type} 1,2,4-triazine ring with nitrogens at positions 1, 2, and 4"

    if not match_found:
        return False, "No valid 1,2,4-triazine ring found"

    return False, "Structure does not match 1,2,4-triazine pattern"