"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem
from rdkit.Chem import rdMolStandardize

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline or its hydrogenated derivative based on its SMILES string.
    A beta-carboline is defined as any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline or hydrogenated derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Beta-carboline SMILES (fully unsaturated)
    beta_carboline_smiles = 'C1=CC=C2C(=C1)C3=C(N2)C=CN=C3'
    beta_carboline_mol = Chem.MolFromSmiles(beta_carboline_smiles)
    if beta_carboline_mol is None:
        return False, "Error generating beta-carboline molecule"

    # Adjust query to ignore bond orders (allowing for hydrogenated derivatives)
    params = rdMolStandardize.AdjustQueryParameters()
    params.makeBondsGeneric = True  # Ignore bond types
    beta_carboline_query = rdMolStandardize.AdjustQueryProperties(beta_carboline_mol, params)

    # Perform substructure search
    if not mol.HasSubstructMatch(beta_carboline_query):
        return False, "Beta-carboline skeleton not found"

    return True, "Contains beta-carboline skeleton or its hydrogenated derivative"