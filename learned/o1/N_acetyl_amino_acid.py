"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    An N-acetyl amino acid is an amino acid where an acetyl group is attached to the nitrogen atom of the amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # SMARTS pattern for acetylated nitrogen
    # Matches any nitrogen atom attached to an acetyl group (N-C(=O)C)
    n_acetyl_pattern = Chem.MolFromSmarts("[N;X3;H0][C;X3](=O)[C;X4H3]")
    acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    if not acetyl_matches:
        return False, "No N-acetyl group found"
    
    # SMARTS pattern for amino acid backbone (generalized)
    # Matches nitrogen connected to alpha-carbon connected to carboxyl group
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[N][C][C](=O)O")
    backbone_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)
    if not backbone_matches:
        return False, "No amino acid backbone found"
    if len(backbone_matches) != 1:
        return False, "More than one amino acid backbone found"

    # Ensure that the acetylated nitrogen is part of the amino acid backbone
    n_idx = backbone_matches[0][0]
    if n_idx not in [match[0] for match in acetyl_matches]:
        return False, "Acetyl group not attached to amino nitrogen"

    return True, "Molecule is an N-acetyl amino acid"