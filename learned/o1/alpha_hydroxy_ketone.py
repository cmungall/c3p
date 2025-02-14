"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-Hydroxy Ketone
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone containing a hydroxy group on 
    the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens (important for checking hydrogen counts)
    mol = Chem.AddHs(mol)

    found_alpha_hydroxy_ketone = False

    # Iterate over atoms to find ketone groups
    for atom in mol.GetAtoms():
        # Check if the atom is a carbonyl carbon (carbon double-bonded to oxygen)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            is_ketone = False
            oxygen_double_bond = None
            carbonyl_neighbors = []
            # Examine bonds of the carbon atom
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Check for double bond to oxygen
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    oxygen_double_bond = neighbor
                # Collect single-bonded neighbors (potential alpha carbons)
                elif bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 6:
                    carbonyl_neighbors.append(neighbor)
            # Confirm it's a ketone (carbonyl carbon bonded to two carbons)
            if oxygen_double_bond and len(carbonyl_neighbors) == 2:
                is_ketone = True

            if is_ketone:
                # Check alpha carbons for hydroxyl groups
                for alpha_carbon in carbonyl_neighbors:
                    # Check if alpha carbon has a hydroxyl group
                    for bond in alpha_carbon.GetBonds():
                        neighbor = bond.GetOtherAtom(alpha_carbon)
                        # Single bond to oxygen atom
                        if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                            oxygen = neighbor
                            # Check if oxygen atom is bonded to a hydrogen (hydroxyl group)
                            num_hydrogens = sum(1 for nbr in oxygen.GetNeighbors() if nbr.GetAtomicNum() == 1)
                            if num_hydrogens > 0:
                                found_alpha_hydroxy_ketone = True
                                return True, "Contains an alpha-hydroxy ketone group"

    if not found_alpha_hydroxy_ketone:
        return False, "Does not contain an alpha-hydroxy ketone group"