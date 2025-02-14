"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide that exhibits antimicrobial properties.
    
    This function checks for features common in peptide antibiotics:
    - Presence of multiple peptide bonds (amide bonds)
    - Molecule is a peptide (contains amino acid residues connected via peptide bonds)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is likely a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define peptide bond pattern (amide bond between C=O and N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    if num_peptide_bonds < 5:
        return False, f"Contains {num_peptide_bonds} peptide bonds; too few for peptide antibiotic"
    
    # Optionally, check for unusual amino acids or modifications
    # For simplicity, we will assume that the presence of multiple peptide bonds is sufficient
    
    return True, f"Contains {num_peptide_bonds} peptide bonds; likely a peptide antibiotic"