"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA basic structure
    adenine = Chem.MolFromSmarts("n1c(nc2c(N)ncnc12)")
    ribose = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)")
    if not (mol.HasSubstructMatch(adenine) and mol.HasSubstructMatch(ribose)):
        return False, "Missing CoA core structure (adenine/ribose)"

    # Check for 4 negative charges (phosphates)
    phosphate = Chem.MolFromSmarts("[O-]P([O-])")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate))
    if phosphate_matches < 2:  # Need at least 2 doubly charged phosphates
        return False, f"Insufficient negative charges, found {phosphate_matches} phosphate groups"

    # Check for 3-oxo group and thioester linkage
    oxo_thioester = Chem.MolFromSmarts("C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_thioester):
        return False, "Missing 3-oxo group with thioester linkage"

    # Check for fatty acid chain (at least 4 carbons)
    fatty_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "Missing sufficient fatty acid chain"

    # Check for pantetheine part
    pantetheine = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[CH](O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Count total carbons in fatty acid chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:  # Rough minimum for CoA + short fatty acid
        return False, f"Carbon count too low ({carbon_count}), expecting at least 25"

    return True, "Contains 3-oxo fatty acid chain connected to CoA(4-) via thioester bond"