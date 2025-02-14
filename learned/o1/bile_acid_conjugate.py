"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a hydrophilic group such as glycine,
    taurine, amino acids, sulfate, glucuronic acid, glucose, other sugars, or coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for bile acid core
    # Steroid nucleus: four fused rings (three 6-membered rings and one 5-membered ring)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC4)C")
    
    # Hydroxyl groups at positions typical for bile acids
    hydroxyl_pattern = Chem.MolFromSmarts("[C;R][C;R](O)[C;R]")
    
    # Check for steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Check for bile acid hydroxylation pattern (optional, can be more specific)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient hydroxyl groups for bile acid"

    # Define SMARTS patterns for conjugated groups
    conjugates = {
        "glycine": Chem.MolFromSmarts("NCC(O)=O"),
        "taurine": Chem.MolFromSmarts("NCCS(=O)(=O)O"),
        "amino_acid": Chem.MolFromSmarts("N[C@@H](C)C(=O)O"),  # Simplified amino acid pattern
        "sulfate": Chem.MolFromSmarts("OS(=O)(=O)[O-]"),
        "glucuronic_acid": Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),
        "glucose": Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
        "coenzyme_A": Chem.MolFromSmarts("NC(=O)CCNC(=O)CC(=O)NCCS"), # Simplified CoA pattern
    }
    
    # Check for conjugation with any of the specified groups
    conjugated = False
    for name, pattern in conjugates.items():
        if mol.HasSubstructMatch(pattern):
            conjugated = True
            conjugate_name = name
            break
    
    if not conjugated:
        return False, "No conjugated hydrophilic group found"

    return True, f"Contains bile acid core conjugated with {conjugate_name}"