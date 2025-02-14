"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom
    which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize a flag to indicate detection of primary alcohol
    found_primary_alcohol = False

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is oxygen
        if atom.GetAtomicNum() == 8:
            # Check if the oxygen has exactly one hydrogen attached (OH group)
            if atom.GetTotalNumHs(includeNeighbors=True) == 1:
                # Iterate over neighbors to oxygen
                for neighbor in atom.GetNeighbors():
                    # Check if neighbor is carbon
                    if neighbor.GetAtomicNum() == 6:
                        carbon = neighbor
                        # Check if carbon is saturated (sp3 hybridized)
                        if carbon.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                            # Count the number of hydrogens attached to carbon
                            num_C_H = carbon.GetTotalNumHs()
                            if num_C_H == 2 or num_C_H == 3:
                                # Count the number of non-hydrogen neighbors excluding oxygen
                                non_H_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() != 1 and nbr.GetIdx() != atom.GetIdx()]
                                num_non_H_neighbors = len(non_H_neighbors)
                                # Check if the carbon is connected to at most one other carbon atom
                                num_C_neighbors = sum(1 for nbr in non_H_neighbors if nbr.GetAtomicNum() == 6)
                                if num_C_neighbors <= 1:
                                    found_primary_alcohol = True
                                    return True, "Contains a primary alcohol group"

    if not found_primary_alcohol:
        return False, "No primary alcohol group found"