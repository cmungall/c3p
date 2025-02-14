"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_alkaloid(smiles: str):
    """
    Determine if a molecule is likely an alkaloid based on its SMILES string.
    This function checks for basic nitrogen in a complex polycyclic structure typical of alkaloids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define alkaloid-relevant nitrogen-containing ring systems
    alkaloid_substructures = [
        Chem.MolFromSmarts('n1ccccc1'),  # Imidazole or pyrrole-like
        Chem.MolFromSmarts('C1=CC=NC=C1'),  # Pyridine, quinoline
        Chem.MolFromSmarts('c1cc[nH]c1'),  # Indole
        Chem.MolFromSmarts('C1CCNCC1')    # Piperidine, pyrrolidine
    ]

    # Check for nitrogen in recognized alkaloid substructures
    for idx, substructure in enumerate(alkaloid_substructures):
        if mol.HasSubstructMatch(substructure):
            return True, f"Contains recognized alkaloid-like nitrogen ring system #{idx + 1}"

    # Additional filtering to avoid common confounders (peptides, specific exocyclic patterns)
    peptide_pattern = Chem.MolFromSmarts('[NX3][CX3](=O)[#6]')
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide-like structure, likely not an alkaloid"

    # Include logic to refine based on other substructures if needed
    # ...

    return False, "No recognizable alkaloid structure found"

# Example usage:
# smiles = "CCCCCCCCCCC1=CC=C(N1)\C=C1/N=C(C=C1OC)C1=CC=CN1"  # Example SMILES
# result, reason = is_alkaloid(smiles)
# print(result, reason)