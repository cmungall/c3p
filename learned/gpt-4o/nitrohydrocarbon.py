"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule contains only C, H, N, O atoms
    # Verify that nitro groups are attached to carbon atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {6, 1, 7, 8}:  # Carbon, Hydrogen, Nitrogen, Oxygen
            return False, f"Contains atom other than C, H, N, O: {atom.GetSymbol()}"

    # Define the nitro group pattern
    nitro_group = Chem.MolFromSmarts("[NX3](=O)[O-]")

    # Find nitro groups in the molecule
    for match in mol.GetSubstructMatches(nitro_group):
        # Check if the first atom in the match (nitrogen in nitro) is connected to a carbon atom
        nitro_nitrogen = mol.GetAtomWithIdx(match[0])
        if any(neighbor.GetAtomicNum() == 6 for neighbor in nitro_nitrogen.GetNeighbors()):
            return True, "Contains hydrocarbon structure with one or more nitro groups replacing hydrogen"
    
    return False, "No nitro group replacing hydrogen found"