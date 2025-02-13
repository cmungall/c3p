"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is defined as a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate net charge of the molecule
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    # Check if net charge is positive
    if net_charge > 0:
        return True, f"Net positive charge of {net_charge} detected"

    return False, "No net positive charge found"