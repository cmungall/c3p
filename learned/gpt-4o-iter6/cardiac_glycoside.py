"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are characterized by a steroid core, lactone ring, and sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for steroid nucleus (cyclopenta[a]phenanthrene)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C1C4CC(C2)CCC4C3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Improved pattern for a 5-membered lactone ring (butenolide)
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Improved pattern for sugar moieties (pyranose and furanose)
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(C(C1O)O)CO"),  # Pyranose
        Chem.MolFromSmarts("C1OC(CO1)CO")       # Furanose
    ]
    sugar_found = any(mol.HasSubstructMatch(sugar) for sugar in sugar_patterns)
    if not sugar_found:
        return False, "No sugar moieties found"

    return True, "Contains steroid nucleus with lactone ring and sugar moieties indicative of cardiac glycoside"