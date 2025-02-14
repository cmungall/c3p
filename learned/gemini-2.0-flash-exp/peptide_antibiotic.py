"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by multiple amide bonds, amino acid structures and a larger MW.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for multiple amide bonds
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:  # Require at least 2 amide bonds
         return False, f"Found {len(amide_matches)} amide bonds, need at least 2"
    
    # Check for alpha carbon
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX2][CX4](C(=O)[OX1])[#6]") # Carbon with N and C=O and another carbon
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) < 2: # Must have at least 2 aminoacids
        return False, f"Found {len(alpha_carbon_matches)} alpha carbon patterns, need at least 2"
    
    # Check for a reasonable number of rotatable bonds to see if it can be a peptide
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5 :
         return False, "Too few rotatable bonds for peptide antibiotic"
    
    # Check for molecular weight (peptide antibiotics are usually > 500)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for peptide antibiotic"
    
    return True, "Contains multiple amide bonds and amino acid structure with a relatively high number of rotatable bonds"