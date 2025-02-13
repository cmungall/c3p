"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:33709 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    Proteinogenic amino acids are the 23 alpha-amino acids that are precursors to proteins
    and are incorporated into proteins during translation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino acid pattern (N-C-C-C=O)
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H1&!$(NC=[!#6])][C@H]([C,c])([C,c])[C](=O)[O,N]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Not an alpha-amino acid"

    # Check for chirality (excluding glycine)
    if mol.GetAtomWithIdx(1).GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        if not mol.GetSmiles() == "NCC(=O)O":  # Glycine
            return False, "Not L-configured"

    # Check for common side chains
    side_chains = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "C(=O)N", "CC(C)C", "CC(N)=O",
                   "CS", "CC1=CC=CC=C1", "CC1=CNC=N1", "CC1=CN=CN1", "C(C)(C)CC(C)(N)C(=O)O"]

    side_chain = Chem.DeleteSubstructs(mol, amino_acid_pattern)
    side_chain_smiles = Chem.MolToSmiles(side_chain)

    if side_chain_smiles not in side_chains:
        if side_chain_smiles == "CC(N)CC(=O)NC(C1=CN=CNC1=O)C(=O)O":  # Pyrrolysine
            return True, "Proteinogenic amino acid (Pyrrolysine)"
        elif side_chain_smiles == "CSCCC":  # Selenocysteine
            return True, "Proteinogenic amino acid (Selenocysteine)"
        else:
            return False, "Unknown side chain"

    return True, "Proteinogenic amino acid"