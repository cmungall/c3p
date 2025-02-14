"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES to create RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino group (-NH2 or protonated forms)
    amino_pattern = Chem.MolFromSmarts("N")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Check for carboxyl group, considering protonated and deprotonated forms (-C(=O)O or -C(=O)[O-])
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"
    
    # Check if glycine pattern matches: it's an exception due to it being non-chiral
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected, which is an achiral proteinogenic amino acid"
    
    # Check for an alpha carbon with chiral configuration connecting both groups
    alpha_carbon_pattern = Chem.MolFromSmarts("C[C@H](N)C(=O)O") 
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "No alpha carbon configuration found with required stereochemistry"
    
    return True, "The molecule fits the profile of a proteinogenic amino acid"