"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33909 corrinoid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond
    linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify pyrrole rings (5-membered rings with one nitrogen)
    pyrrole_rings = []
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    for ring in atom_rings:
        if len(ring) == 5:
            num_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if num_nitrogens == 1:
                pyrrole_rings.append(ring)
    
    if len(pyrrole_rings) < 4:
        return False, f"Found {len(pyrrole_rings)} pyrrole rings, need at least 4"
    
    # Attempt to determine if the pyrrole rings are connected in a macrocycle
    # This is complex to determine accurately
    # We can attempt to check if there is a large ring containing the pyrrole rings
    
    # Find all ring systems (connected rings)
    # First, mark the atoms in the pyrrole rings
    pyrrole_atoms = set()
    for ring in pyrrole_rings:
        pyrrole_atoms.update(ring)
    
    # Check for a macrocyclic ring involving these atoms
    sssr = Chem.GetSymmSSSR(mol)
    macrocycle_found = False
    for ring in sssr:
        if len(ring) >= 15:  # Assuming macrocycle is at least 15 atoms
            if pyrrole_atoms.issubset(set(ring)):
                macrocycle_found = True
                break
    
    if not macrocycle_found:
        return False, "No macrocyclic ring containing the pyrrole rings found"
    
    # Check for cobalt atom (many corrinoids contain cobalt coordinated to the nitrogen atoms)
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if not has_cobalt:
        return False, "No cobalt atom found in the molecule"
    
    # At this point, we can classify the molecule as a corrinoid
    return True, "Molecule contains corrin nucleus with four pyrrole rings in a macrocycle"