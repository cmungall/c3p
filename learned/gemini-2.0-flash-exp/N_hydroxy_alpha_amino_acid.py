"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    An N-hydroxy-alpha-amino acid is an amino acid in which at least one hydrogen 
    attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for the alpha amino acid core: a carboxyl group and an alpha carbon bonded to a nitrogen
    alpha_amino_acid_core = Chem.MolFromSmarts("[NX3]-[CX4]-[CX3](=[OX1])-[OX2]")
    if not mol.HasSubstructMatch(alpha_amino_acid_core):
        return False, "No alpha amino acid core found"

    # Pattern for N-hydroxy group attached to the amino nitrogen: -N(OH)-C-C(=O)O or  -N(OH)(OH)-C-C(=O)O
    n_hydroxy_pattern_1 = Chem.MolFromSmarts("[NX3]([OX2])-[CX4]-[CX3](=[OX1])-[OX2]")  # -N(OH)-C-C(=O)O
    n_hydroxy_pattern_2 = Chem.MolFromSmarts("[NX3]([OX2])([OX2])-[CX4]-[CX3](=[OX1])-[OX2]")  # -N(OH)(OH)-C-C(=O)O
    n_hydroxy_pattern_3 = Chem.MolFromSmarts("[NX2]=[CX3]-[OX2]") # -N=C-OH amidoxime
    
    if mol.HasSubstructMatch(n_hydroxy_pattern_1) or mol.HasSubstructMatch(n_hydroxy_pattern_2) or mol.HasSubstructMatch(n_hydroxy_pattern_3):
        return True, "Contains an alpha amino acid structure with at least one hydroxy group on the amino nitrogen"
    
    return False, "No N-hydroxy group attached to the alpha amino group found"