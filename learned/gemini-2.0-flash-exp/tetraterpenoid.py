"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are C40 isoprenoids, often with modifications like methyl loss/rearrangement.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check carbon count. Must be around 40.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 38 or c_count > 44:
        return False, f"Carbon count {c_count} is not in range. Must be between 38 and 44"

    # Check presence of isoprene units (C5 unit: CC(=C)C).
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 5:
        return False, f"Too few isoprene units: {len(isoprene_matches)}. Requires at least 5"

    # Check presence of double bonds. Tetraterpenoids are highly unsaturated.
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3 :
       return False, f"Too few double bonds: {len(double_bond_matches)}. Requires at least 3"


    # Check molecular weight. Tetraterpenoids tend to be heavy > 500. Allow flexibility.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
       return False, f"Molecular weight too low: {mol_wt}. Requires at least 500"

    return True, "Meets tetraterpenoid criteria: C40 isoprenoid structure, multiple double bonds."