"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core oxime substructure using SMARTS (C=N-O)
    oxime_pattern = Chem.MolFromSmarts("[C]=[N][O]")
    
    # Find matches for the oxime pattern
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
      return False, "Molecule does not contain the C=N-O oxime group"

    # Check if at least one of the matched C is [CH1]
    for match in matches:
        c_atom_idx = match[0] #The carbon in the C=N
        c_atom = mol.GetAtomWithIdx(c_atom_idx)
        if c_atom.GetTotalNumHs() == 1: #if this C is a [CH1]
            return True, "Molecule contains the R-CH=NOH aldoxime substructure."

    return False, "Molecule does not contain the correct R-CH=NOH aldoxime substructure."