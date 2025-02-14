"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide - A peptide containing ten or more amino acid residues.
"""

from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide containing ten or more amino acid residues)
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the peptide bond SMARTS pattern (peptide bonds are amide bonds between amino acids)
    peptide_bond_smarts = "C(=O)N"

    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)

    if peptide_bond_pattern is None:
        return False, "Invalid peptide bond SMARTS pattern"

    # Find peptide bonds in the molecule
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)

    # Number of amino acid residues is number of peptide bonds + 1
    num_amino_acids = num_peptide_bonds + 1

    if num_amino_acids >= 10:
        return True, f"Contains {num_amino_acids} amino acid residues connected via peptide bonds"
    else:
        return False, f"Only {num_amino_acids} amino acid residues, less than 10 required for polypeptide"