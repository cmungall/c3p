"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids,
    typically between 2 and 20 residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond SMARTS pattern (amide bond between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond_pattern is None:
        return False, "Invalid SMARTS pattern for peptide bond"

    # Find peptide bonds (amide bonds)
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    # Count the number of peptide bonds
    num_peptide_bonds = len(peptide_bonds)

    # Estimate the number of amino acid residues
    # For linear peptides, number of residues = number of peptide bonds + 1
    # For cyclic peptides, this estimation may differ
    num_amino_acids = num_peptide_bonds + 1 if num_peptide_bonds > 0 else 0

    # Check if there are any peptide bonds
    if num_amino_acids == 0:
        return False, "No peptide bonds found"

    # Check if the number of residues is within oligopeptide range (2-20)
    if num_amino_acids < 2:
        return False, f"Molecule has {num_amino_acids} amino acid residue(s), fewer than required for an oligopeptide"
    if num_amino_acids > 20:
        return False, f"Molecule has {num_amino_acids} amino acid residues, more than allowed for an oligopeptide"

    # Additional check to exclude simple amide-containing molecules
    # Verify that the molecule has amino acid-like substructures
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # Simplified amino acid pattern
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid residues detected"

    return True, f"Molecule has {num_amino_acids} amino acid residues, classified as an oligopeptide"