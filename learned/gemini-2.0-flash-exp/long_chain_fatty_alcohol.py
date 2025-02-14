"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13-C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty alcohol, False otherwise.
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of alcohol groups (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    num_hydroxyls = len(hydroxyl_matches)
    if num_hydroxyls != 1:
        return False, f"Molecule has {num_hydroxyls} alcohol groups, it should have exactly 1"


    # Count the number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check if the chain length is within the range of C13-C22
    if carbon_count < 13 or carbon_count > 22:
          return False, f"Molecule has {carbon_count} carbons, should be between 13 and 22"

    # Check that the only heteroatoms are a single oxygen from the -OH
    heteroatoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1]
    if len(heteroatoms) != 1 or heteroatoms[0] != 8: # Only one heteroatom and must be Oxygen
       return False, "Molecule contains heteroatoms other than a single oxygen"

    return True, "Molecule is a long-chain fatty alcohol (C13-C22) with one alcohol group."