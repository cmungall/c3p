"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol is a compound in which a hydroxy group (-OH) is attached to
    a saturated carbon atom which has two other carbon atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found_secondary_alcohol = False

    # Iterate over all atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            if atom.GetDegree() == 1:  # Connected to only one atom
                neighbor = atom.GetNeighbors()[0]
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    # Check if the carbon is sp3 hybridized (saturated)
                    if neighbor.GetHybridization() == rdchem.HybridizationType.SP3:
                        # Count the number of carbon neighbors
                        carbon_neighbors = 0
                        for nbr in neighbor.GetNeighbors():
                            if nbr.GetAtomicNum() == 6:
                                carbon_neighbors += 1
                        if carbon_neighbors == 2:
                            found_secondary_alcohol = True
                            break  # Found a secondary alcohol group

    if found_secondary_alcohol:
        return True, "Contains a secondary alcohol functional group"
    else:
        return False, "No secondary alcohol functional group found"