"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:15554 dipeptide

A dipeptide is any molecule that contains two amino-acid residues connected by peptide linkages.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bond pattern (-C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Count number of peptide bonds
    n_peptide_bonds = len(peptide_bond_matches)
    if n_peptide_bonds != 2:
        return False, f"Found {n_peptide_bonds} peptide bonds, expected 2 for a dipeptide"

    # Look for amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H1][C;H1]([C;H2,H1,H0])[C;H2,H1,H0]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Count number of amino acid residues
    n_amino_acids = len(amino_acid_matches)
    if n_amino_acids != 2:
        return False, f"Found {n_amino_acids} amino acid residues, expected 2 for a dipeptide"

    # Check for terminal groups (-COOH and -NH2)
    carboxyl_pattern = Chem.MolFromSmarts("[C;H1](=O)[O;H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    amine_pattern = Chem.MolFromSmarts("[N;H2]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    if len(carboxyl_matches) != 1 or len(amine_matches) != 1:
        return False, "Missing terminal -COOH or -NH2 group"

    return True, "Contains two amino acid residues connected by peptide linkages"