"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a benzene ring with exactly four chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
         return False, "No benzene ring found"


    # Find all chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    
    # count how many chlorine atoms are directly attached to the benzene
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) > 1:
         return False, "More than one benzene ring found"
    
    benzene_atoms_indices = benzene_matches[0]
    
    chlorine_on_benzene_count = 0
    for cl_atom in chlorine_atoms:
        for neighbor_atom in cl_atom.GetNeighbors():
            if neighbor_atom.GetIdx() in benzene_atoms_indices:
                 chlorine_on_benzene_count += 1

    if chlorine_on_benzene_count != 4:
        return False, f"Found {chlorine_on_benzene_count} chlorine substituents on the benzene ring, need exactly 4"

    return True, "Benzene ring with exactly four chlorine substituents found"