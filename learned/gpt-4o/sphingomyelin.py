"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelin is characterized by a sphingoid base linked to a fatty acid and phosphorylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key structural features of sphingomyelin
    amide_bond_pattern = Chem.MolFromSmarts("N-C(=O)-C")
    phosphorylcholine_pattern = Chem.MolFromSmarts("COP(OCC[N+](C)(C)C)(=O)[O-]")

    # Check for amide bond
    if not mol.HasSubstructMatch(amide_bond_pattern):
        return False, "No amide bond with fatty acid found"

    # Check for phosphorylcholine
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "No phosphorylcholine headgroup found"

    return True, "Contains amide-linked fatty acid and phosphorylcholine headgroup"