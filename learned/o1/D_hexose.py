"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16551 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a monosaccharide with six carbons and D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Find all carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) != 6:
        return False, f"Number of carbon atoms is {len(carbons)}, expected 6 for hexose"

    # Find aldehyde group (C=O) to identify open-chain form
    aldehyde_smarts = '[CX3H1](=O)'
    aldehyde = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde)
    if aldehyde_matches:
        # Open-chain form detected
        c1_idx = aldehyde_matches[0][0]
        # Traverse the chain to number the carbons from C1 to C6
        carbon_chain = [c1_idx]
        visited = set(carbon_chain)
        current_idx = c1_idx

        while len(carbon_chain) < 6:
            atom = mol.GetAtomWithIdx(current_idx)
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                         if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
            if neighbors:
                current_idx = neighbors[0]
                carbon_chain.append(current_idx)
                visited.add(current_idx)
            else:
                break  # No further carbons found

        if len(carbon_chain) != 6:
            return False, "Could not find a chain of 6 connected carbons starting from the aldehyde carbon"

        # Get the C5 atom (5th carbon in the chain)
        c5_idx = carbon_chain[4]
        c5_atom = mol.GetAtomWithIdx(c5_idx)

        # Check if C5 is chiral and get its configuration
        if not c5_atom.HasProp('_CIPCode'):
            return False, "C5 is not a chiral center"
        c5_chirality = c5_atom.GetProp('_CIPCode')
        if c5_chirality == 'R':
            return True, "Molecule is a D-hexose; C5 has R configuration"
        elif c5_chirality == 'S':
            return False, "Molecule is an L-hexose; C5 has S configuration"
        else:
            return False, "Chirality at C5 is undefined"

    else:
        # Cyclic form detected; classification not handled
        return None, "Cyclic forms are not handled by this classifier"