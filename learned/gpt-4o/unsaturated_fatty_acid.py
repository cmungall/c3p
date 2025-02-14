"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a terminal carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of terminal carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No terminal carboxyl group found (-COOH)"
    
    # Check for presence of unsaturation: C=C bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond_pattern):
        unsat_bonds = True
    else:
        # Check for presence of unsaturation: C#C bonds
        triple_bond_pattern = Chem.MolFromSmarts("C#C")
        unsat_bonds = mol.HasSubstructMatch(triple_bond_pattern)
    
    if not unsat_bonds:
        return False, "No unsaturation (C=C or C#C) found in the structure"

    # Check for long aliphatic chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:  # Assuming a minimum chain length of 10 carbons
        return False, "Carbon chain too short to be considered a fatty acid"

    return True, "Contains at least one unsaturation and terminal carboxyl group, typical of unsaturated fatty acids"