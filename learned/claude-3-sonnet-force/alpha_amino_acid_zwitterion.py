"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:49804 - Alpha-amino-acid zwitterion
An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; 
major species at pH 7.3.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for zwitterionic pattern: NH3+ and COO-
    zwitterion_pattern = Chem.MolFromSmarts("[NH3+][CX3](C([O-])=O)[H]")
    
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No zwitterionic NH3+/COO- group found"
    
    # Check for alpha position (alpha carbon bound to both NH3+ and COO-)
    alpha_carbon_atoms = mol.GetSubstructMatches(zwitterion_pattern)
    alpha_carbon_atom = alpha_carbon_atoms[0][1]
    
    # Count number of substituents on alpha carbon (excl. NH3+ and COO-)
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_atom)
    substituents_on_alpha = sum(1 for atom in alpha_carbon.GetNeighbors() if atom.GetSymbol() != "N" and atom.GetSymbol() != "O")
    
    if substituents_on_alpha != 1:
        return False, "Alpha carbon is not singly substituted (not an alpha-amino acid)"
    
    return True, "Contains alpha-amino-acid zwitterionic structure"