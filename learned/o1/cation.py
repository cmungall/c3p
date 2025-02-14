"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton,
    resulting in a net positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cation, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate the net formal charge of the molecule
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    if net_charge > 0:
        return True, f"Molecule has net positive charge ({net_charge}); it is a cation"
    else:
        return False, f"Molecule has net charge ({net_charge}); it is not a cation"