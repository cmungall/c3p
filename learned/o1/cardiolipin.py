"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin consists of two phosphatidic acid units connected via a central glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count phosphorus atoms (P)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 2:
        return False, f"Expected 2 phosphorus atoms, found {p_count}"

    # Count ester bonds (C(=O)O-C)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_ester_bonds = len(ester_matches)
    if num_ester_bonds < 4:
        return False, f"Expected at least 4 ester bonds (fatty acid chains), found {num_ester_bonds}"

    # Check for central glycerol unit connected via phosphodiester bonds
    # Central glycerol connected to two phosphate groups
    central_glycerol_pattern = Chem.MolFromSmarts('[C;!R]-[C;!R]-[C;!R]')
    glycerol_matches = mol.GetSubstructMatches(central_glycerol_pattern)
    if not glycerol_matches or len(glycerol_matches) < 1:
        return False, "Central glycerol unit not found"

    # Check for phosphodiester bonds connected to central glycerol
    phosphodiester_pattern = Chem.MolFromSmarts('P(=O)(O)(O)-O-C')
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if len(phosphodiester_matches) < 2:
        return False, f"Expected at least 2 phosphodiester bonds, found {len(phosphodiester_matches)}"

    # Check for phosphatidic acid units (glycerol with two esterified fatty acids and one phosphate group)
    phosphatidic_acid_pattern = Chem.MolFromSmarts('C(CO[P](=O)(O)O)(COC(=O)C)OC(=O)C')
    phosphatidic_acid_matches = mol.GetSubstructMatches(phosphatidic_acid_pattern)
    if len(phosphatidic_acid_matches) < 2:
        return False, f"Expected at least 2 phosphatidic acid units, found {len(phosphatidic_acid_matches)}"

    # Molecular weight check (cardiolipins are large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low for cardiolipin ({mol_wt} Da)"

    # All checks passed
    return True, "Molecule matches cardiolipin structural features"