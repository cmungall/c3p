"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is formed by the esterification of the carboxylic acid group
    of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group pattern (ester linkage)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    if ester_pattern is None:
        return False, "Invalid ester SMARTS pattern"

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group found"

    # Define the alpha-amino acid ester pattern
    # Nitrogen attached to alpha carbon, which is attached to carbonyl carbon of ester
    amino_acid_ester_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CH1,CH2][CX3](=O)[OX2H0][#6]")
    if amino_acid_ester_pattern is None:
        return False, "Invalid amino acid ester SMARTS pattern"

    # Check for the alpha-amino acid ester pattern
    if mol.HasSubstructMatch(amino_acid_ester_pattern):
        return True, "Molecule contains an alpha-amino acid ester moiety"

    # Alternatively, check for alpha carbon connected to amino group and ester
    for ester_match in ester_matches:
        ester_carbon_idx = ester_match[0]  # Index of the carbonyl carbon in the ester
        ester_carbon = mol.GetAtomWithIdx(ester_carbon_idx)

        # Find alpha carbons (connected to ester carbonyl carbon)
        alpha_carbons = [atom for atom in ester_carbon.GetNeighbors() if atom.GetAtomicNum() == 6]
        for alpha_carbon in alpha_carbons:
            # Check if alpha carbon is connected to an amino group
            has_amino_group = False
            for neighbor in alpha_carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 7:
                    # Ensure the nitrogen is an amino nitrogen (not part of an amide or other group)
                    if neighbor.GetDegree() <= 3 and not neighbor.IsInRing():
                        has_amino_group = True
                        break
            if has_amino_group:
                return True, "Molecule contains an alpha-amino acid ester moiety"

    return False, "Molecule does not contain an alpha-amino acid ester moiety"