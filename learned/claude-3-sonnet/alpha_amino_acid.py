"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:16669 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is an amino acid in which the amino group is located on the carbon atom
    at the position alpha to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;H2][C;X4][C;X3](=O)[O;H1]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        # Assign stereochemistry and check if the matched substructure is alpha
        AllChem.AssignStereochemistry(mol)
        matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
        for match in matches:
            alpha_carbon = mol.GetAtomWithIdx(match[1])
            if alpha_carbon.GetHybridization() == Chem.HybridizationType.SP3:
                return True, "Contains an amino group at the alpha position to a carboxyl group"

    return False, "No alpha-amino acid substructure found"