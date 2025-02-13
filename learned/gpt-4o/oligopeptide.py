"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is usually a short chain of amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an oligopeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to create an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for peptide bonds considering possible alpha-carbon and nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3]")
    if peptide_bond_pattern is None:
        return False, "Unable to create SMARTS pattern for peptide bond"

    # Find matches for the peptide bond pattern
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # Check if the number of peptide bonds is consistent with an oligopeptide
    if 2 <= num_peptide_bonds < 20:
        return True, f"Contains {num_peptide_bonds} peptide bonds, consistent with an oligopeptide"
    elif num_peptide_bonds >= 20:
        return False, f"Contains {num_peptide_bonds} peptide bonds, too many for an oligopeptide"
    else:
        return False, "No peptide bonds or insufficient for classification as an oligopeptide"