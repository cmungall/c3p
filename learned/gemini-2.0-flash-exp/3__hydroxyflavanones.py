"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone core with a hydroxy substituent at the 3' position
    of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flavanone core pattern. X means any atom, [CH] is a specific atom.
    flavanone_core_pattern = Chem.MolFromSmarts("C1=CC=C2[C](C(C=O)CC2)[CH]1") # modified to look for the chiral center
    if not mol.HasSubstructMatch(flavanone_core_pattern):
          return False, "Molecule does not contain flavanone core"


    # 3'-hydroxy substitution pattern - this searches the 3' position of the B ring. The ring is numbered clockwise starting with the carbon attached to the flavanone.
    hydroxy_3_pattern = Chem.MolFromSmarts("c1[c]([OH])[c](c[c]c1)-[CH]") # the carbon attached to the ring should be linked to the core structure

    if not mol.HasSubstructMatch(hydroxy_3_pattern):
        return False, "Molecule does not contain a hydroxyl at 3' position"


    return True, "Molecule contains a flavanone core with a hydroxyl group at the 3' position"