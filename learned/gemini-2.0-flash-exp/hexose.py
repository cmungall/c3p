"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
       return False, f"Molecule does not have exactly 6 carbons, it has {c_count}."

    # Check for pyranose (6-membered ring with O)
    pyranose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4][CX4]1")
    pyranose_matches = mol.HasSubstructMatch(pyranose_pattern)

    # Check for furanose (5-membered ring with O)
    furanose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4]1")
    furanose_matches = mol.HasSubstructMatch(furanose_pattern)


    # If it's not a ring, then check for linear chain with a carbonyl
    if not pyranose_matches and not furanose_matches:
        # Check for aldehyde at C1
        aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)[C]") # C(=O) next to another carbon
        aldehyde_matches = mol.HasSubstructMatch(aldehyde_pattern)

         #Check for ketone at C2
        ketone_pattern = Chem.MolFromSmarts("[CX4][C](=O)[CX4]") # ketone (C=O) connected to two carbons.
        ketone_matches = mol.HasSubstructMatch(ketone_pattern)

        if not aldehyde_matches and not ketone_matches:
            return False, "Molecule is not a ring or a linear chain with a carbonyl at the 1 or 2 position."


    return True, "Molecule is a hexose (6 carbons with aldehyde at C1 or ketone at C2 or ring structure)"