"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: CHEBI:36164 Lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy substituents (-OO).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for lipid backbone (long carbon chain)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 12:
        return False, "Not enough carbons for a lipid backbone"

    # Look for hydroperoxy groups (-OO)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy groups found"

    # Check for carboxylic acid group (lipid terminus)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found (lipid terminus)"

    # Check for unsaturated bonds (typical for lipid hydroperoxides)
    n_unsaturated = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetBondType() == Chem.BondType.TRIPLE)
    if n_unsaturated < 2:
        return False, "Too few unsaturated bonds for a lipid hydroperoxide"

    # Count rotatable bonds to verify lipid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Lipid chain too short or rigid"

    return True, "Contains lipid backbone with one or more hydroperoxy groups"