"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:36354 tannin

A tannin is defined as any of a group of astringent polyphenolic vegetable principles or compounds,
chiefly complex glucosides of catechol and pyrogallol. Key features:
- Polyphenolic compounds
- Often contain glucose or gallic acid derivatives
- Multiple phenol/catechol/pyrogallol rings connected by ester or ether linkages
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for polyphenolic structure (multiple phenol/catechol/pyrogallol rings)
    phenol_pattern = Chem.MolFromSmarts("[c;H1]")  # Aromatic ring with 1 attached H
    catechol_pattern = Chem.MolFromSmarts("[c;H0;r5]1:c:c:c:c:c:1")  # Aromatic ring with 2 attached O
    pyrogallol_pattern = Chem.MolFromSmarts("[c;H0;r5]1:c:c(:O):c:c:c:1")  # Aromatic ring with 3 attached O

    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pattern)

    if len(phenol_matches) + len(catechol_matches) + len(pyrogallol_matches) < 3:
        return False, "Fewer than 3 phenol/catechol/pyrogallol rings found"

    # Look for glucose or gallic acid derivatives
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H]([C@@H]([C@@H]1O)O)O)CO")
    gallic_pattern = Chem.MolFromSmarts("C(=O)c1cc(O)c(O)c(O)c1")

    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    gallic_matches = mol.GetSubstructMatches(gallic_pattern)

    if not glucose_matches and not gallic_matches:
        return False, "No glucose or gallic acid derivatives found"

    # Look for ester/ether linkages between rings
    ester_pattern = Chem.MolFromSmarts("[O;X2]-[CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("[O;X2]-[#6]")

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    if not ester_matches and not ether_matches:
        return False, "No ester or ether linkages found between rings"

    return True, "Polyphenolic compound containing glucose/gallic acid derivatives and phenol/catechol/pyrogallol rings connected by ester/ether linkages"