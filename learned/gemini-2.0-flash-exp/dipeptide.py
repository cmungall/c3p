"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule consisting of two amino acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Dipeptide pattern: two amino acids connected by one peptide bond.
    dipeptide_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])N[CX4][CX3](=[OX1])[OX2]")
    
    if not mol.HasSubstructMatch(dipeptide_pattern):
        return False, "Molecule does not contain the dipeptide substructure (two amino acids connected by peptide bond)"

    matches = mol.GetSubstructMatches(dipeptide_pattern)
    if len(matches) != 1:
        return False, f"Molecule has {len(matches)} dipeptide substructures, should have exactly 1."

    return True, "Contains two amino acid residues connected by a single peptide bond"