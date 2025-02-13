"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic linkages, commonly involving sugar units like beta-D-glucopyranosides
    # We'll define a simple SMARTS pattern for a generic glucopyranoside
    glucopyranoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](CO)O1")
    if not mol.HasSubstructMatch(glucopyranoside_pattern):
        return False, "No sugar moiety (glucoside) found"

    # Check for triterpenoid structure, which would typically involve a C30 terpene structure
    # With multiple rings and functional groups, such as hydroxyl or acetyl groups.
    # Example SMARTS pattern for lupane/oleanane type triterpenes (very simplified)
    triterpenoid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC5C4CC(CC(C5)C)(C)O")
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "No triterpenoid moiety found"

    # Check for linkage between glycoside and triterpenoid
    # Often saponins are linked through an ester or ether bond
    linkage_pattern = Chem.MolFromSmarts("O[*]-C=O")  # Simplified ester pattern
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No linkage between glycoside and triterpenoid"

    # Ensure the overall structure has a significant molecular weight (typically >800 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for typical triterpenoid saponin"

    return True, "Contains characteristics of a triterpenoid saponin with glycoside linkage"