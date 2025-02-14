"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI:49686 lipopeptide
A compound consisting of a peptide with attached lipid.
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import FragmentMatcher

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide contains both a peptide and a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for peptide backbone (N-C=O)
    peptide_pattern = Chem.MolFromSmarts("[N;X3][C;X3](=[O;X1])")
    has_peptide = mol.HasSubstructMatch(peptide_pattern)
    
    # Look for lipid chains (long aliphatic groups)
    lipid_pattern = Chem.MolFromSmarts("[C;X4][C;X4]~[C;X4][C;X4]~[C;X4][C;X4]~[C;X4][C;X4]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    has_lipid = len(lipid_matches) > 0
    
    # Classify as lipopeptide if both components are present
    if has_peptide and has_lipid:
        return True, "Contains peptide backbone and lipid chains"
    elif has_peptide:
        return False, "Contains peptide backbone but no lipid chains"
    elif has_lipid:
        return False, "Contains lipid chains but no peptide backbone"
    else:
        return False, "No peptide or lipid components found"