"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are long-chain polyunsaturated fatty acids with multiple
    cis double bonds and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for carboxylic acid group (COOH) pattern at one end of the molecule
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for multiple double bonds (C=C) along the chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bonds) < 3:
        return False, "Too few double bonds; polyunsaturation required"

    # Verify the presence of cis configuration double bonds
    cis_bond_count = mol.GetSubstructMatches(Chem.MolFromSmarts('C/C=C/C')) + mol.GetSubstructMatches(Chem.MolFromSmarts('C\C=C\C'))
    if len(cis_bond_count) < len(double_bonds):
        return False, "Not enough cis double bonds for essential fatty acid"

    # Verify the total number of carbon atoms in the chain
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18:
        return False, "Too few carbons for essential fatty acid"

    return True, "Contains carboxylic acid and sufficient cis double bonds for essential fatty acid"