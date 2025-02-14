"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for 1-O-acylglycerophosphoethanolamine components
    
    # 1-O-acylglycerol pattern: includes stereochemistry for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO[*])[CH2]O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1-O-acyl-glycerol structure found"

    # Ester bond pattern at 1-position (attached to the glycerol)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[C@@H]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond at 1-position of glycerol found"

    # Phosphoethanolamine group pattern
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not present"

    # Verify oxygen positions; one from ester, two from phosphoester
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 7:  # The necessary minimum number of oxygens
        return False, f"Insufficient oxygen atoms, found {oxygen_count}"

    return True, "Structure matches 1-O-acylglycerophosphoethanolamine"