"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol with a single nitro group attached directly on the same benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible SMARTS for phenol with a single nitro group on benzene ring 
    phenol_nitro_pattern = Chem.MolFromSmarts("c1(ccccc1O)[NX3](=O)[O-]")
    
    # Look for a phenol ring with a single nitro group
    matches = mol.GetSubstructMatches(phenol_nitro_pattern)
    if len(matches) == 0:
        return False, "No phenol with a single nitro group found on the benzene ring"
    
    # Check that there is exactly one nitro group on the ring
    nitro_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N' and any(neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors())]
    if len(nitro_groups) != 1:
        return False, f"Invalid number of nitro groups: found {len(nitro_groups)}"
    
    return True, "Contains phenol group with a single nitro group on the same benzene ring"