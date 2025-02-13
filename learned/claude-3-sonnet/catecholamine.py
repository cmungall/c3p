"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18237 catecholamine
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a 4-(2-aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core pyrocatechol (benzene-1,2-diol) structure
    pyrocatechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)ccc1")
    if not mol.HasSubstructMatch(pyrocatechol_pattern):
        return False, "No pyrocatechol moiety found"

    # Look for the 2-aminoethyl side chain (-CH2CH2NH2)
    aminoethyl_pattern = Chem.MolFromSmarts("NCCCC")
    aminoethyl_match = mol.GetSubstructMatches(aminoethyl_pattern)
    if not aminoethyl_match:
        return False, "No 2-aminoethyl side chain found"

    # Check if the side chain is attached to the pyrocatechol ring
    for match in aminoethyl_match:
        atom_idx = match[-1]  # index of the carbon atom attached to the ring
        if mol.GetAtomWithIdx(atom_idx).IsInRingSize(6):
            return True, "Contains the core 4-(2-aminoethyl)pyrocatechol structure"

    return False, "2-aminoethyl side chain not attached to pyrocatechol ring"