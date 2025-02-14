"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find primary amine nitrogen atoms (degree 1 nitrogen with no charge)
    amine_atoms = [atom for atom in mol.GetAtoms()
                   if atom.GetAtomicNum() == 7 and atom.GetDegree() == 1 and atom.GetFormalCharge() == 0]
    if not amine_atoms:
        return False, "No primary amino group found"

    for amine_atom in amine_atoms:
        n_idx = amine_atom.GetIdx()

        # Get the carbon atom attached to the amine nitrogen (C2)
        c2_atoms = [nbr for nbr in amine_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not c2_atoms:
            continue  # No carbon attached to nitrogen
        c2_atom = c2_atoms[0]
        c2_idx = c2_atom.GetIdx()

        # Check for hydroxyl groups on carbons adjacent to C2 (C1 and C3)
        has_adjacent_hydroxyl = False
        
        # Get neighbors of C2 atom
        c2_neighbors = c2_atom.GetNeighbors()
        for atom in c2_neighbors:
            if atom.GetAtomicNum() == 6 and atom.GetIdx() != n_idx:  # Exclude nitrogen
                # Check if this carbon has a hydroxyl group
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                        has_adjacent_hydroxyl = True
                        break
            if has_adjacent_hydroxyl:
                break

        if not has_adjacent_hydroxyl:
            continue  # No hydroxyl group adjacent to amino group

        # Now, check the aliphatic chain length starting from C2
        visited = set()
        chain_atoms = []

        def traverse_chain(atom):
            if atom.GetIdx() in visited:
                return
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() != 6:
                return
            chain_atoms.append(atom)
            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.IsInRing():
                    continue  # Skip ring structures
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != amine_atom.GetIdx():
                    traverse_chain(nbr)

        traverse_chain(c2_atom)
        chain_length = len(chain_atoms)

        if chain_length < 12:
            continue  # Chain too short to be sphingoid

        # If all checks passed
        return True, "Molecule matches sphingoid structural features"

    return False, "Molecule does not match sphingoid structural features"