"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has an amino group and a carboxylate group attached to the same
    chiral alpha carbon, which is also connected to another substituent (e.g., alkyl, aryl, hydroxyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino-acid zwitterion pattern
    zwitterion_pattern = Chem.MolFromSmarts("[C@H]([NH3+])([C,c])(C(=O)[O-])")
    zwitterion_matches = mol.GetSubstructMatches(zwitterion_pattern)
    if not zwitterion_matches:
        return False, "No alpha-amino-acid zwitterion pattern found"

    # Check for overall charge neutrality
    mol_charge = AllChem.GetFormalCharge(mol)
    if mol_charge != 0:
        return False, "Molecule is not charge-neutral"

    # Additional checks for alpha-amino-acid zwitterions
    # ...

    return True, "Molecule is an alpha-amino-acid zwitterion"