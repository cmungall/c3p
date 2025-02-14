"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid in which the amino acid specified is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define N-acylglycine pattern
    # The pattern represents R-C(=O)-N-CH2-C(=O)-O
    n_acylglycine_pattern = Chem.MolFromSmarts('O=CNCC(=O)O')
    if not mol.HasSubstructMatch(n_acylglycine_pattern):
        return False, "No N-acylglycine substructure found"

    # Ensure that the nitrogen is part of an amide linkage (acylated)
    amide_pattern = Chem.MolFromSmarts('C(=O)NCC(=O)O')
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Nitrogen is not acylated"
    
    # Check that there are no substituents on the CH2 group (glycine backbone)
    glycine_pattern = Chem.MolFromSmarts('NCC(=O)O')
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    if not glycine_matches:
        return False, "No glycine backbone found"
    for match in glycine_matches:
        # Index of the CH2 carbon
        ch2_idx = match[1]
        ch2_atom = mol.GetAtomWithIdx(ch2_idx)
        # Check that the CH2 carbon is only bonded to nitrogen and carbonyl carbon
        neighbors = ch2_atom.GetNeighbors()
        if len(neighbors) != 2:
            return False, "Alpha carbon has substituents, not glycine"
        for nbr in neighbors:
            if nbr.GetAtomicNum() != 7 and nbr.GetAtomicNum() != 6:
                return False, "Alpha carbon bonded to atoms other than N and C"
    
    return True, "Molecule is an N-acylglycine"