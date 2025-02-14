"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate includes a polyprenol chain attached via a phosphate ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate or diphosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    diphosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OP(=O)(O)O")
    if not (mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(diphosphate_pattern)):
        return False, "No phosphate or diphosphate group found"
    
    # Check for ester linkage to phosphate group
    ester_linkage_pattern = Chem.MolFromSmarts("O-P(=O)(O)O")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "Ester linkage to phosphate group not found"

    # Look for repeating isoprene units (C=C-C-C-C pattern)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C")
    isoprene_count = len(mol.GetSubstructMatches(isoprene_pattern))
    if isoprene_count < 3:  # We need multiple units to classify as polyprenol
        return False, f"Insufficient isoprene units: {isoprene_count}"

    return True, "Contains polyprenol chain with phosphate ester bond"