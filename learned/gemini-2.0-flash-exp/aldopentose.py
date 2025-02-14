"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a five-carbon monosaccharide with an aldehyde group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Molecule has {c_count} carbons, but aldopentose has 5."

    # 2. Check for aldehyde group (-C(=O)H or its cyclic equivalent -C(OH)-O-)
    # Open-chain aldehyde (H-C=O)
    aldehyde_open_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])")
    # Cyclic hemiacetal formed with aldehyde
    aldehyde_closed_pattern = Chem.MolFromSmarts("[CX4H1][OX2][CX4]")
    
    if not mol.HasSubstructMatch(aldehyde_open_pattern) and not mol.HasSubstructMatch(aldehyde_closed_pattern):
        return False, "No aldehyde group or cyclic hemiacetal found."

    # 3. Check for 5 oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 5:
      return False, f"Molecule has {o_count} oxygens, but aldopentose has 5."

    #If the molecule passes all the tests then it is classified as an aldopentose
    return True, "Molecule is classified as an aldopentose"