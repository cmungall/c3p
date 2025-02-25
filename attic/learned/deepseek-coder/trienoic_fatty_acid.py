"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:36586 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with exactly three double bonds in the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count double bonds in the carbon chain (excluding carbonyl and ring double bonds)
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Exclude double bonds in rings
            if bond.GetBeginAtom().IsInRing() or bond.GetEndAtom().IsInRing():
                continue
            # Exclude carbonyl double bonds (C=O)
            if bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8:
                continue
            double_bond_count += 1

    if double_bond_count != 3:
        return False, f"Found {double_bond_count} double bonds in the carbon chain, need exactly 3"

    # Check if the molecule is a fatty acid (long carbon chain)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, "Carbon chain too short to be a fatty acid"

    return True, "Contains a carboxylic acid group and exactly three double bonds in the carbon chain"