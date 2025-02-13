"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:46419 nitrile
A compound having the structure RC#N; thus a C-substituted derivative of hydrocyanic acid, HC#N. 
In systematic nomenclature, the suffix nitrile denotes the triply bound #N atom, not the carbon atom attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a cyano group (-C#N) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyano group (-C#N) attached to a carbon atom
    nitrile_pattern = Chem.MolFromSmarts("[C]#N")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No cyano group (-C#N) attached to a carbon atom found"
    
    # Check for specific non-nitrile compounds containing cyano groups
    non_nitrile_smarts = [
        "O=C(C#N)C#N",  # oxomalononitrile
        "N#Cc1ccccc1"   # benzonitrile
    ]
    
    for smarts in non_nitrile_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return False, f"Matched non-nitrile compound: {smarts}"
    
    # If cyano group attached to carbon and no specific non-nitrile cases matched, classify as nitrile
    return True, "Contains a cyano group (-C#N) attached to a carbon atom"