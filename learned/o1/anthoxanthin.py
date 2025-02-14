"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments in plants, including flavones and flavonols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavone core SMILES (2-phenylchromen-4-one)
    flavone_core_smiles = 'O=C1Oc2ccccc2C=C1c3ccccc3'
    flavone_core_mol = Chem.MolFromSmiles(flavone_core_smiles)
    if flavone_core_mol is None:
        return False, "Failed to create flavone core molecule"

    # Define flavonol core SMILES (3-hydroxyflavone)
    flavonol_core_smiles = 'O=C1Oc2ccc(O)cc2C=C1c3ccccc3'
    flavonol_core_mol = Chem.MolFromSmiles(flavonol_core_smiles)
    if flavonol_core_mol is None:
        return False, "Failed to create flavonol core molecule"

    # Check for substructure match with flavone core
    if mol.HasSubstructMatch(flavone_core_mol):
        return True, "Molecule contains flavone core structure"

    # Check for substructure match with flavonol core
    if mol.HasSubstructMatch(flavonol_core_mol):
        return True, "Molecule contains flavonol core structure"

    # If no match is found
    return False, "Molecule does not contain flavone or flavonol core structure"