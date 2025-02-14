"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:32988 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PeptideDrawing

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

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

    # Check if molecule is a peptide
    peptide = PeptideDrawing.DetectPeptide(mol)
    if not peptide:
        return False, "Not a peptide"

    # Count number of amino acid residues
    residues = PeptideDrawing.DetectPeptideResidues(mol)
    n_residues = len(residues)

    # Oligopeptides typically have < 10 amino acids
    if n_residues < 10:
        return True, f"Contains {n_residues} amino acid residues"
    else:
        return False, f"Contains {n_residues} amino acid residues, likely not an oligopeptide"