"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: Polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is a member of the class of phenols that contains two or more benzene rings,
    each of which is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure proper aromaticity perception
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        # Handle molecules that RDKit cannot kekulize
        Chem.Kekulize(mol, clearAromaticFlags=True)

    # Exclude peptides by checking for peptide bonds (amide linkages)
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains peptide bonds, not a polyphenol"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Set to keep track of benzene rings with hydroxy substitution
    benzene_rings_with_hydroxy = set()

    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Not a six-membered ring

        # Check if all atoms in the ring are carbons
        is_all_carbon = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)
        if not is_all_carbon:
            continue  # Ring contains non-carbon atoms

        # Check if all atoms in the ring are aromatic
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_aromatic:
            continue  # Ring is not aromatic

        # Check if any atom in the ring has a hydroxy substitution
        has_hydroxy = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                nbr_atom = neighbor
                # Check for hydroxy group (Oxygen atom connected to hydrogen and the ring carbon)
                if nbr_atom.GetAtomicNum() == 8 and nbr_atom.GetDegree() == 1:
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break  # No need to check other atoms in this ring

        if has_hydroxy:
            benzene_rings_with_hydroxy.add(frozenset(ring))

    ring_count = len(benzene_rings_with_hydroxy)

    if ring_count >= 2:
        return True, f"Molecule is a polyphenol with {ring_count} benzene rings each substituted by at least one hydroxy group"
    else:
        return False, f"Molecule has {ring_count} benzene rings with hydroxy substitution, needs at least 2 to be a polyphenol"