"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorinated biphenyl (PCB) based on its SMILES string.
    A PCB is a biphenyl compound containing between 2 and 10 chlorine atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PCB, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for biphenyl structure (two connected benzene rings)
    biphenyl_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[c]2[c][c][c][c][c]2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Not a biphenyl structure"
    
    # Count chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)

    # Check if chlorine count is within the range 2-10
    if chlorine_count < 2:
         return False, f"Too few chlorine atoms: {chlorine_count}. Must be between 2 and 10."
    if chlorine_count > 10:
        return False, f"Too many chlorine atoms: {chlorine_count}. Must be between 2 and 10."


    return True, "Biphenyl structure with 2-10 chlorine atoms"