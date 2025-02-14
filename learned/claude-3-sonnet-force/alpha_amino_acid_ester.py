"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:51366 alpha-amino acid ester
Alpha-amino acid esters are the amino acid ester derivatives obtained by the formal condensation of an alpha-amino acid with an alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find amino acid backbone (-C(=O)N-) with at least one attached carbon
    aa_pattern = Chem.MolFromSmarts("[C](=[O])([N])([C])")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    
    if len(aa_matches) == 0:
        return False, "No alpha-amino acid backbone found"
    
    # Check for ester group (-C(=O)O-) attached to the alpha carbon
    ester_pattern = Chem.MolFromSmarts("[C]([C](=[O])([O]))")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) == 0:
        return False, "No ester group found on alpha carbon"
    
    # Check if ester carbon is part of amino acid backbone
    for ester_match in ester_matches:
        if ester_match[0] in [atom_idx for match in aa_matches for atom_idx in match]:
            return True, "Contains alpha-amino acid backbone with an ester group on the alpha carbon"
    
    return False, "Ester group not part of alpha-amino acid backbone"