"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to match L-alpha-amino acid core structure with L-configuration
    #   [#6] : Carbon atom
    #   [C@@H] : Carbon atom with L-stereoconfiguration
    #   [NX3;H2,H3;!+] : Nitrogen with three bonds (either -NH2 or -NH3+)
    #   [CX3](=[OX1]) : Carboxyl group (-COOH or -COO-)
    #   ~ : Indicates any kind of bond
    #   [!#6] any atom that is not C (for the carboxyl)
    #   [!#1] any atom that is not hydrogen (for the amino group)

    # SMARTS for the L-alpha-amino acid core (neutral and zwitterionic forms)
    smarts_neutral = "[NX3;H2;!+][C@@H]([!#1])([CX3](=[OX1])[!#6])"
    smarts_zwitter = "[NX3;H3;!+][C@@H]([!#1])([CX3](=[OX1])[OX1-])"
    smarts_zwitter2 = "[NX3;H3;!+][C@@H]([!#1])([CX3](=O)[O-])"
    # Check if the molecule matches the SMARTS pattern

    core_pattern_neutral = Chem.MolFromSmarts(smarts_neutral)
    core_pattern_zwitter = Chem.MolFromSmarts(smarts_zwitter)
    core_pattern_zwitter2 = Chem.MolFromSmarts(smarts_zwitter2)

    if not (mol.HasSubstructMatch(core_pattern_neutral) or mol.HasSubstructMatch(core_pattern_zwitter) or mol.HasSubstructMatch(core_pattern_zwitter2)):
         return False, "Molecule does not contain the L-alpha-amino acid core structure."
    

    return True, "Molecule contains the L-alpha-amino acid core structure with L-configuration."