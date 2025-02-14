"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule that contains two amino acid residues connected by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (amide bond between two amino acids)
    peptide_bond_smarts = '[NX3][CX3](=O)'
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)

    # Find all peptide bonds
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond)
    num_peptide_bonds = len(peptide_bond_matches)

    if num_peptide_bonds != 1:
        return False, f"Found {num_peptide_bonds} peptide bonds, need exactly 1 for dipeptide"

    # Define amino acid residue pattern (alpha amino acid backbone)
    # N-C(alpha)-C(=O)
    amino_acid_smarts = '[NX3][CX4H1][CX3](=O)'
    amino_acid = Chem.MolFromSmarts(amino_acid_smarts)

    # Find all amino acid residues
    amino_acid_matches = mol.GetSubstructMatches(amino_acid)
    num_amino_acids = len(amino_acid_matches)

    if num_amino_acids != 2:
        return False, f"Found {num_amino_acids} amino acid residues, need exactly 2"

    # Check if the amino acids are connected via the peptide bond
    # Map the atoms involved in peptide bond
    peptide_bond_atoms = [match for match in peptide_bond_matches[0]]
    amino_acid_atoms = [list(match) for match in amino_acid_matches]

    # Verify that the peptide bond connects the two amino acid residues
    connected = False
    for aa1 in amino_acid_atoms:
        for aa2 in amino_acid_atoms:
            if aa1 != aa2:
                # Check if peptide bond connects aa1 and aa2
                if (peptide_bond_atoms[0] in aa1 and peptide_bond_atoms[1] in aa2) or \
                   (peptide_bond_atoms[0] in aa2 and peptide_bond_atoms[1] in aa1):
                    connected = True
                    break
        if connected:
            break

    if not connected:
        return False, "Peptide bond does not connect the two amino acid residues"

    return True, "Contains two amino acid residues connected by a peptide bond"