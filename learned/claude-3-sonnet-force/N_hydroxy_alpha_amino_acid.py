"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: CHEBI:26205 N-hydroxy-alpha-amino-acid
Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino acid backbone pattern (C-C-N-C-C-C=O)
    amino_acid_pattern = Chem.MolFromSmarts("[C,$(C(-C)(-C)(-C))](-C(-C(-N(-[OH])))(-C(-C(-C(=O)))))")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"
    
    # Look for N-hydroxy group (-N(-[OH])-)
    n_hydroxy_pattern = Chem.MolFromSmarts("N(-[OH])")
    if not mol.HasSubstructMatch(n_hydroxy_pattern):
        return False, "No N-hydroxy group found"
    
    # Check for alpha position of N-hydroxy group
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetDegree() == 3:
            n_hydroxy_atom = atom
            break
    else:
        return False, "N-hydroxy group not in alpha position"
    
    if n_hydroxy_atom.GetTotalNumHs() != 1:
        return False, "N-hydroxy group not in alpha position"
    
    return True, "Contains N-hydroxy group in alpha position of amino acid backbone"