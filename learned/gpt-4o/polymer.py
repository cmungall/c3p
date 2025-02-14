"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Focuses on identifying complex chains and repeating units, common in polymers.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is identified as a polymer-like structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect typical polymer-like repeating units or diverse elements
    patterns = [
        Chem.MolFromSmarts("[C](=[O,N])-[C]"),  # common in polyesters, polyamides
        Chem.MolFromSmarts("[C]([C])[C]([C])[C]"),  # common in polyethylene chains
        Chem.MolFromSmarts("[O]-[C]-[C]-[O]"),  # polyether linkage
    ]

    # Check matches for polymer-like chains
    total_matches = 0
    for pattern in patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_matches += len(matches)

    if total_matches > 10:  # Aiming for a larger occurrence indicating polymer
        return True, f"Contains repeating substructural patterns: found {total_matches}"

    # Analyze molecular size indicative of polymer
    if mol.GetNumAtoms() > 100:  # refined larger size threshold
        return True, f"Large molecular size with {mol.GetNumAtoms()} atoms, possible polymer"

    # Count number of single bond chains greater than a certain length indicative of backbone
    if max(rdmolops.GetDistanceMatrix(mol).flatten()) > 15:  # Longest path indicative of polymer chains
        return True, "Contains long chains indicative of polymeric structure"

    return False, "No explicit polymer-like features detected"