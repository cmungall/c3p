"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:24869 cation
A monoatomic or polyatomic species having one or more elementary charges of the proton.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for positive formal charges on atoms
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge > 0:
            return True, f"Contains atom with positive formal charge: {formal_charge}"
    
    # Look for bracketed cationic elements (e.g. [Na+], [NH4+])
    cationic_pattern = Chem.MolFromSmarts("[+]")
    if mol.HasSubstructMatch(cationic_pattern):
        return True, "Contains bracketed cationic element"
    
    # Look for common organic cation substructures
    cation_patterns = [
        Chem.MolFromSmarts("[NH3+]"),  # Ammonium
        Chem.MolFromSmarts("[N+]"),    # Quaternary ammonium
        Chem.MolFromSmarts("[S+]"),    # Sulfonium
        Chem.MolFromSmarts("[P+]"),    # Phosphonium
        Chem.MolFromSmarts("[O+]"),    # Oxonium
    ]
    for pattern in cation_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains organic cation substructure: {Chem.MolToSmiles(pattern)}"
    
    return False, "No evidence of cationic character"