"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:36915 cation
A monoatomic or polyatomic species having one or more elementary charges of the proton.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species with one or more positive charges.

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

    # Calculate formal charge on the molecule
    formal_charge = Chem.GetFormalCharge(mol)

    # Check for positive formal charge
    if formal_charge > 0:
        return True, f"Molecule has a positive formal charge of {formal_charge}"
    else:
        return False, "Molecule does not have a positive formal charge"

    # Alternative approach: Check for positive charges on atoms
    # has_positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    # if has_positive_charge:
    #     return True, "Molecule contains positively charged atoms"
    # else:
    #     return False, "No positively charged atoms found"