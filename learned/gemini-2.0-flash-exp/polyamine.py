"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as an organic compound with two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for different types of amino groups. Explicitly look for N with H to differentiate it from amides.
    amino_patterns = [
        Chem.MolFromSmarts("[NH2;X4]"),  # Primary amines (with explicit hydrogens)
        Chem.MolFromSmarts("[NH;X4;H1]"),   # Secondary amines
        Chem.MolFromSmarts("[N;X4;H0]"),   # Tertiary amines
        Chem.MolFromSmarts("[NH3+]"),         # Primary protonated amines
        Chem.MolFromSmarts("[NH2+]"),          # Secondary protonated amines
        Chem.MolFromSmarts("[nH]1[c,n,o,s]~[c,n,o,s]~[c,n,o,s]1"), #Ring amines,
         Chem.MolFromSmarts("[c,C]~[NH2;X4]"), # Aromatic Amines

        
    ]

    total_amino_count = 0
    matched_atoms = set()

    for pattern in amino_patterns:
        if pattern is None: continue
        matches = mol.GetSubstructMatches(pattern)

        for match in matches:
            for atom_index in match:
                if atom_index not in matched_atoms:
                    total_amino_count += 1
                    matched_atoms.add(atom_index)


    if total_amino_count >= 2:
        return True, "Contains two or more amino groups"
    else:
        return False, "Does not contain two or more amino groups"