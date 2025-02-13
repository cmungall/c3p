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

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for non-hydrocarbon atoms, excluding nitro group atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {6, 1, 7, 8}:  # Carbon, Hydrogen, Nitrogen, Oxygen
            return False, f"Contains atom other than C, H, N, O: {atom.GetSymbol()}"

    # Check if the molecule has at least one nitro group directly attached to a hydrocarbon structure
    nitro_group = Chem.MolFromSmarts("[NX3](=O)[O-]")  # Nitro group pattern
    for match in mol.GetSubstructMatches(nitro_group):
        # Check if Nitro group is attached to a carbon
        if any(mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors()):            
            return True, "Contains hydrocarbon structure with one or more nitro groups replacing hydrogen"
    
    return False, "No nitro group replacing hydrogen found"