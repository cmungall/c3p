"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:51443 anilide

An anilide is any aromatic amide obtained by acylation of aniline.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aniline substructure (aromatic ring with NH2 group)
    aniline_pattern = Chem.MolFromSmarts("c1ccc(cc1)N")
    aniline_matches = mol.GetSubstructMatches(aniline_pattern)

    # Look for amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    # Check if aniline and amide groups are connected
    connected = False
    for aniline_match in aniline_matches:
        for amide_match in amide_matches:
            aniline_atom = aniline_match[0]
            amide_atom = amide_match[1]
            if mol.GetBondBetweenAtoms(aniline_atom, amide_atom):
                connected = True
                break

    if connected:
        return True, "Contains aniline group connected to an amide group"
    else:
        return False, "Aniline and amide groups are not connected"