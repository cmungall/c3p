"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:17164 primary amine
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound where a nitrogen atom is bonded to one hydrocarbyl group and two hydrogens.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Atomic number for nitrogen
            if atom.GetFormalCharge() != 0:
                continue  # Skip if nitrogen has a formal charge
                
            # Total number of hydrogens (implicit + explicit)
            total_hydrogens = atom.GetTotalNumHs()
            
            # Continue if nitrogen has exactly two hydrogens
            if total_hydrogens != 2:
                continue

            # Get neighbors of the nitrogen atom
            neighbors = atom.GetNeighbors()
            
            # Initialize counters
            num_carbon = 0
            num_heteroatoms = 0
            
            for neighbor in neighbors:
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 6:
                    num_carbon += 1
                elif atomic_num != 1:  # Not hydrogen
                    num_heteroatoms += 1

            # Check if nitrogen is bonded to exactly one carbon and no other heteroatoms
            if num_carbon == 1 and num_heteroatoms == 0:
                return True, "Primary amine group found"

    return False, "No primary amine group found"