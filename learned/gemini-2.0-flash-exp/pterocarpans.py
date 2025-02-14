"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan has a 6H-[1]benzofuro[3,2-c]chromene skeleton.
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the core 6H-[1]benzofuro[3,2-c]chromene skeleton with the 6a and 11a positions
    # '[CH1]1[CH2]2Oc3c4ccccc4Oc4c3cc5ccccc5[CH1]21'
    core_pattern = Chem.MolFromSmarts('[CH1]1[CH2]2Oc3c4ccccc4Oc4c3cc5ccccc5[CH1]21')

    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 6H-[1]benzofuro[3,2-c]chromene skeleton not found."

    # Verify 6a and 11a positions are dihydro
    match = mol.GetSubstructMatch(core_pattern)
    if match:
         atom_6a = mol.GetAtomWithIdx(match[0])
         atom_11a = mol.GetAtomWithIdx(match[len(match)-1])
         if (atom_6a.GetTotalNumHs() != 1 or atom_11a.GetTotalNumHs() !=1):
            return False, "6a and 11a are not dihydro"


    return True, "Molecule contains the core 6H-[1]benzofuro[3,2-c]chromene skeleton with dihydro 6a and 11a positions"