"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl has 2-10 chlorine atoms attached to a biphenyl core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for biphenyl core with more flexible pattern
    # This pattern matches two aromatic rings connected by a single bond
    biphenyl_pattern = Chem.MolFromSmarts("c:1:c:c:c:c:c:1-c:1:c:c:c:c:c:1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"

    # Count chlorine atoms
    cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    n_chlorines = len(cl_atoms)
    
    if n_chlorines < 2:
        return False, f"Too few chlorines ({n_chlorines}), minimum is 2"
    if n_chlorines > 10:
        return False, f"Too many chlorines ({n_chlorines}), maximum is 10"

    # Get all atoms in the molecule
    all_atoms = mol.GetAtoms()
    
    # Check for unwanted elements or substituents
    allowed_atoms = {6, 17}  # Only carbon and chlorine allowed
    for atom in all_atoms:
        atomic_num = atom.GetAtomicNum()
        # Skip hydrogen as it's often implicit
        if atomic_num == 1:
            continue
        if atomic_num not in allowed_atoms:
            return False, f"Contains non-allowed element: {atom.GetSymbol()}"
        
        # For carbon atoms, check they're either part of aromatic system or single bonds only
        if atomic_num == 6:
            if not atom.GetIsAromatic():
                # Check all bonds are single
                bonds = atom.GetBonds()
                for bond in bonds:
                    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        return False, "Contains non-aromatic carbon with non-single bonds"

    # Get the two benzene rings from biphenyl core
    matches = mol.GetSubstructMatches(biphenyl_pattern)
    if len(matches) > 0:
        ring_atoms = set(matches[0])
        
        # Check if all atoms are part of the core structure or are chlorines
        for atom in all_atoms:
            idx = atom.GetIdx()
            if idx not in ring_atoms and atom.GetAtomicNum() != 17:
                return False, "Contains atoms outside biphenyl core that aren't chlorine"

        # Verify chlorines are attached to aromatic carbons
        for cl in cl_atoms:
            neighbors = cl.GetNeighbors()
            if len(neighbors) != 1:
                return False, "Invalid chlorine bonding"
            neighbor = neighbors[0]
            if not (neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic()):
                return False, "Chlorine not attached to aromatic carbon"
            if neighbor.GetIdx() not in ring_atoms:
                return False, "Chlorine attached outside biphenyl core"

    return True, f"Valid polychlorobiphenyl with {n_chlorines} chlorine atoms"