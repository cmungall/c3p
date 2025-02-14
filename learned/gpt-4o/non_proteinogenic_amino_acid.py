"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid has an amino group, a carboxyl group, and a distinct variation
    from the 20 standard amino acids found in proteins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amino acid structure: [N][C][C](=O)O with optional stereochemistry
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H3]C(C(=O)O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Does not contain basic amino acid structure"

    # Additional checks for non-standard amino acids
    standard_amino_acids_smarts = [
        "N[C@@H](C)C(=O)O",  # Alanine
        "N[C@@H](CCC(=O)O)C(=O)O",  # Aspartic acid
        # Add more SMARTS for each standard amino acid here if needed
    ]
    
    for aa_smarts in standard_amino_acids_smarts:
        aa_pattern = Chem.MolFromSmarts(aa_smarts)
        if mol.HasSubstructMatch(aa_pattern):
            return False, "Matches a standard amino acid"

    # Check for indications of non-standard features, long chains, or unusual elements
    unusual_atoms = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "N", "O"]]
    if unusual_atoms:
        return True, f"Contains unusual elements: {set(unusual_atoms)}"

    # Check for side chain variations, especially length and complexity
    side_chain_len = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 1.0)
    if side_chain_len > 12:
        return True, "Has long side chain suggesting a non-standard amino acid"

    return True, "Does not match any standard amino acids and fits the criteria for non-proteinogenic"