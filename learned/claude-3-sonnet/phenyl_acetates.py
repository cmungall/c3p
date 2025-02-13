"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
           with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a phenyl acetate and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for phenyl acetate:
    # More specific pattern requiring phenolic oxygen and acetate group
    phenyl_acetate_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1][cR1][cR1][cR1]1-[OX2][CX3](=[OX1])[CH3]")
    
    # Negative patterns - structures that should not be considered phenyl acetates
    negative_patterns = [
        # Complex ring systems
        Chem.MolFromSmarts("[R2]"), # Any atom in 2 or more rings
        # Specific exclusions
        Chem.MolFromSmarts("O=C1Oc2ccccc2C1"),  # coumarin
        Chem.MolFromSmarts("C1=NCc2ccccc12"),    # indole
        Chem.MolFromSmarts("O=C1C=Cc2ccccc2O1"), # benzofuranone
        Chem.MolFromSmarts("O=C1C=COc2ccccc21"), # benzopyranone
        # Exclude certain functional groups that indicate more complex structures
        Chem.MolFromSmarts("[N+](=O)[O-]"),      # nitro group
        Chem.MolFromSmarts("[#7R]"),             # nitrogen in ring
        Chem.MolFromSmarts("[OR2]"),             # oxygen in multiple rings
    ]

    # First check negative patterns
    for pattern in negative_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains excluded structural features"

    # Check for phenyl acetate pattern
    matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if not matches:
        return False, "No phenyl acetate group found"

    # Additional validation for each match
    for match in matches:
        # Get the phenyl ring atoms
        ring_atoms = match[:6]
        
        # Verify ring is a single benzene ring (not part of larger system)
        ring_info = mol.GetRingInfo()
        ring_count = 0
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count rings this atom is part of
            ring_count += sum(1 for ring in ring_info.AtomRings() if atom_idx in ring)
            
        # Each atom should be in exactly one ring
        if ring_count != 6:
            continue

        # Verify acetate group
        oxygen_idx = match[6]
        carbonyl_idx = match[7]
        methyl_idx = match[9]
        
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        methyl = mol.GetAtomWithIdx(methyl_idx)
        
        # Verify proper connectivity
        if (len(oxygen.GetNeighbors()) != 2 or
            len(carbonyl.GetNeighbors()) != 3 or
            len(methyl.GetNeighbors()) != 1):
            continue
            
        # If we get here, we have a valid match
        return True, "Contains phenyl ring with acetate group properly attached to phenolic oxygen"

    return False, "No valid phenyl acetate group found"