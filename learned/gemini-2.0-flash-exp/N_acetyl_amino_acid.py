"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amine group is acetylated
    (CH3-C(=O)-N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for N-acetyl amino acid with more flexible matching of N and carboxylic acid
    acetyl_amino_acid_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3][CX4][CX3](=[OX1])[OX1,OX2]")

    if mol.HasSubstructMatch(acetyl_amino_acid_pattern):
        # Further check: there should be one N-acetyl and one amino acid core
        acetyl_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3]")
        amino_acid_core = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[OX1,OX2]")
        
        acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
        amino_matches  = mol.GetSubstructMatches(amino_acid_core)


        if len(acetyl_matches) == 1 and len(amino_matches) == 1:
            return True, "Molecule matches the N-acetyl-amino acid pattern."
        else:
            return False, f"Molecule has {len(acetyl_matches)} acetyl group and {len(amino_matches)} amino acid core."
    else:
        return False, "Molecule does not match the N-acetyl-amino acid pattern."