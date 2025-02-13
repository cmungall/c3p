"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group linked to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-acetyl group pattern
    acetyl_pattern = Chem.MolFromSmarts("N[C,CX3](=O)C")  # Allow broader context of acetyl linkage
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No N-acetyl group found"

    # Look for amino acid backbone or variants: allow some flexibility on amino group to catch more configurations
    amino_acid_patterns = [
        Chem.MolFromSmarts("[NX3][C@@H](C)C(=O)O"),   # Common chiral structures
        Chem.MolFromSmarts("[NX3][C@H](C)C(=O)O"),    # Non-chiral counterpart
        Chem.MolFromSmarts("N[*]CC(=O)O"),            # Generic structural variant for non-chiral centers
        Chem.MolFromSmarts("N1[C@@H](C)C1C(=O)O"),    # Cyclic structures like proline
    ]
    
    has_amino_acid_structure = any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns)
    if not has_amino_acid_structure:
        return False, "No identifiable amino acid structure"

    # Check if there's only a single N-acetyl amino acid like structure
    # This checks for overall molecule size constraints to avoid complex structures like peptides
    if mol.GetNumAtoms() > 25:  # Arbitrary limit, might adjust after further tests
        return False, "Molecule too large to be a single N-acetyl-amino acid"

    return True, "Contains N-acetyl group attached to an amino acid structure"