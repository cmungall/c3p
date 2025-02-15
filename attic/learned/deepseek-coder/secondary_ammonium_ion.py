"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:137982 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is an organic cation obtained by protonation of any secondary amino compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure pattern for a secondary ammonium ion (R2NH2+)
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NX3;H2;+1]([#6])([#6])")
    
    # Check if the molecule contains the secondary ammonium ion pattern
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        return True, "Contains a secondary ammonium ion (R2NH2+)"
    else:
        return False, "No secondary ammonium ion (R2NH2+) found"