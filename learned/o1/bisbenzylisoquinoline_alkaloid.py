"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges,
    direct carbon-carbon bridges, or methylenedioxy groups.

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

    # Define benzylisoquinoline unit as SMARTS pattern
    # This pattern represents an isoquinoline ring attached to a benzyl group
    benzylisoquinoline_smarts = "c1ccccc1CC2=CC=NC3=CC=CC=C32"
    benzylisoquinoline_pattern = Chem.MolFromSmarts(benzylisoquinoline_smarts)

    if benzylisoquinoline_pattern is None:
        return False, "Invalid benzylisoquinoline SMARTS pattern"

    # Find all benzylisoquinoline units in the molecule
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} benzylisoquinoline units, need at least 2"

    # Convert matches into list of sets for easier processing
    unit_atom_sets = [set(match) for match in matches]

    # Check for bridges between benzylisoquinoline units
    bridges_found = False
    for i in range(len(unit_atom_sets)):
        for j in range(i+1, len(unit_atom_sets)):
            unit_i_atoms = unit_atom_sets[i]
            unit_j_atoms = unit_atom_sets[j]
            # Iterate over bonds to find connections between units
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                # Check if bond connects atoms from different units
                if (atom1_idx in unit_i_atoms and atom2_idx in unit_j_atoms) or \
                   (atom1_idx in unit_j_atoms and atom2_idx in unit_i_atoms):
                    atom1 = mol.GetAtomWithIdx(atom1_idx)
                    atom2 = mol.GetAtomWithIdx(atom2_idx)
                    # Check if the bridge involves oxygen or carbon atoms
                    if atom1.GetAtomicNum() in [6, 8] and atom2.GetAtomicNum() in [6, 8]:
                        # Check for methylenedioxy group (O-CH2-O)
                        if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6:
                            for nbr in atom2.GetNeighbors():
                                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [atom1_idx, atom2_idx]:
                                    bridges_found = True
                                    return True, "Contains two benzylisoquinoline units linked by methylenedioxy bridge"
                        else:
                            bridges_found = True
                            return True, "Contains two benzylisoquinoline units linked by bridges"
    if bridges_found:
        return True, "Contains two benzylisoquinoline units linked by bridges"
    else:
        return False, "No bridges found between benzylisoquinoline units"