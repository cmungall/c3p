"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: iridoid monoterpenoid
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically consists of a cyclopentane ring fused to a six-membered
    oxygen-containing heterocycle (tetrahydropyran ring), forming a specific bicyclic system.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) == 0:
        return False, "No rings found in molecule"

    # Find all five-membered rings containing only carbons
    five_membered_rings = []
    for ring in ssr:
        if len(ring) == 5:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_nums = [atom.GetAtomicNum() for atom in atoms_in_ring]
            if all(num == 6 for num in atom_nums):
                five_membered_rings.append(set(ring))

    if not five_membered_rings:
        return False, "No cyclopentane ring found"

    # Find all six-membered rings containing one oxygen atom
    six_membered_rings = []
    for ring in ssr:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_nums = [atom.GetAtomicNum() for atom in atoms_in_ring]
            if atom_nums.count(8) == 1 and atom_nums.count(6) == 5:
                six_membered_rings.append(set(ring))

    if not six_membered_rings:
        return False, "No tetrahydropyran ring found"

    # Check for fused rings (shared atoms)
    found_iridoid_core = False
    for ring5 in five_membered_rings:
        for ring6 in six_membered_rings:
            shared_atoms = ring5 & ring6
            if len(shared_atoms) >= 2:
                # Rings are fused
                found_iridoid_core = True
                break
        if found_iridoid_core:
            break

    if not found_iridoid_core:
        return False, "No fused cyclopentane and tetrahydropyran rings found"

    # Verify monoterpenoid backbone (typically at least 10 carbon atoms)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 10:
        return False, f"Too few carbons for monoterpenoid (found {num_carbons}, need at least 10)"

    return True, "Contains iridoid core structure typical of iridoid monoterpenoids"