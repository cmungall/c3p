"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # D-ribose component pattern (matches ribose moiety in nucleosides)
    ribose_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O")
    
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "D-ribose component not found"

    # Check for an attached nucleobase pattern (purine or pyrimidine) - indicative but generic
    purine_base_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")  # Generic purine-like structure
    pyrimidine_base_pattern = Chem.MolFromSmarts("n1c[nH]c(=O)nc1")  # Generic pyrimidine-like structure
    
    has_nucleobase = mol.HasSubstructMatch(purine_base_pattern) or mol.HasSubstructMatch(pyrimidine_base_pattern)

    if not has_nucleobase:
        return False, "Nucleobase component not discernible"

    return True, "Molecule contains a D-ribose moiety linked to a nucleobase, classified as a ribonucleoside"