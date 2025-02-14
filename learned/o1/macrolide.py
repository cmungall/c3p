"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring structures found in molecule"

    # Define SMARTS pattern for lactone group (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C(=O)O[C;R]")
    if lactone_pattern is None:
        return None, "Error in defining lactone SMARTS pattern"

    # Initialize flag for macrolide detection
    is_macrolide_flag = False

    # Iterate over each ring to check for macrocyclic lactone
    for ring_atoms in atom_rings:
        ring_size = len(ring_atoms)
        if ring_size >= 12:
            # Create a substructure of the ring
            ring_mol = Chem.PathToSubmol(mol, ring_atoms)
            # Check if the ring contains a lactone group
            if ring_mol.HasSubstructMatch(lactone_pattern):
                is_macrolide_flag = True
                break  # Found a macrocyclic lactone ring

    if is_macrolide_flag:
        return True, "Contains macrocyclic lactone ring of 12 or more members"
    else:
        return False, "No macrocyclic lactone rings of required size found"