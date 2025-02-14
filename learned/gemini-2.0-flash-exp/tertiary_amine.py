"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all N atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    if not nitrogen_atoms:
      return False, "No nitrogen atoms found"
    
    # check if each N is bonded to exactly three C's
    found_tertiary_amine = False
    for atom in nitrogen_atoms:
        carbon_neighbors = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
            elif neighbor.GetAtomicNum() == 7 or neighbor.GetAtomicNum() == 8 or neighbor.GetAtomicNum() == 16 or neighbor.GetAtomicNum() == 1: # N,O,S or H
                pass
            else:
                 return False, "Nitrogen atom is not only attached to carbon." # Added check for other heteroatoms/H
        
        if carbon_neighbors == 3:
            found_tertiary_amine = True
        elif carbon_neighbors < 3 or carbon_neighbors > 3:
            return False, "Nitrogen atom is not a tertiary amine"

    if found_tertiary_amine:
       return True, "Tertiary amine"
    else:
        return False, "No tertiary amine found"