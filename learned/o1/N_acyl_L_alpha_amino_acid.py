"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is any L-alpha-amino acid carrying an N-acyl substituent.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generalized SMARTS pattern for N-acyl-L-alpha-amino acid
    # [C@@H]: L-stereochemistry at alpha carbon
    # N[C](=O)[#6]: N-acyl group where acyl can be any carbon chain
    # [C](=O)O: Carboxylic acid group
    pattern = Chem.MolFromSmarts('[C@@H](N[C](=O)[#6])[*]C(=O)O')

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        # Check for additional amide bonds to exclude peptides
        amide_pattern = Chem.MolFromSmarts('C(=O)N')
        amide_bonds = mol.GetSubstructMatches(amide_pattern)
        num_amide_bonds = len(amide_bonds)

        if num_amide_bonds > 1:
            return False, f"Found {num_amide_bonds} amide bonds, molecule may be a peptide"
        
        return True, "Molecule is an N-acyl-L-alpha-amino acid"

    # If not found, check without stereochemistry
    pattern_no_stereo = Chem.MolFromSmarts('C(N[C](=O)[#6])[*]C(=O)O')
    if mol.HasSubstructMatch(pattern_no_stereo):
        # Check for additional amide bonds
        amide_pattern = Chem.MolFromSmarts('C(=O)N')
        amide_bonds = mol.GetSubstructMatches(amide_pattern)
        num_amide_bonds = len(amide_bonds)

        if num_amide_bonds > 1:
            return False, f"Found {num_amide_bonds} amide bonds, molecule may be a peptide"
        
        return True, "Molecule is an N-acyl-alpha-amino acid (without specific stereochemistry)"

    return False, "Molecule is not an N-acyl-L-alpha-amino acid"