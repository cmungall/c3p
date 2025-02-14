"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is a peptide consisting of three amino acid residues connected by two peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an alpha carbon in a peptide bond
    # This looks for a carbon with 2 non-hydrogen neighbors. 
    # It must be connected to a nitrogen and also to a carbonyl carbon
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=[OX1])")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    # Count the number of alpha carbons and verify that there are three
    if len(alpha_carbon_matches) != 3:
        return False, f"Found {len(alpha_carbon_matches)} alpha carbons, expected 3"
    
    # check that we have two peptide bonds; we look for a carbon with 3 neighbors, one is oxygen
    # and the other is a nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("[NX2,NX3][CX3](=[OX1])")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
         return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected 2"
    
    return True, "Contains three amino acid residues connected by two peptide bonds."