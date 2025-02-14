"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    This class of molecules has a glycerol backbone with:
    - An acyl chain attached via ester bond at position 1
    - A free hydroxyl group at position 2
    - A phosphoethanolamine group at position 3

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns without chiral specifications
    # Phosphoethanolamine group (allowing for variable protonation states)
    phosphoethanolamine_pattern = Chem.MolFromSmarts('O[P](=O)(O)OCCN')
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # Glycerol backbone with acyl group at position 1 via ester bond
    acyl_glycerol_pattern = Chem.MolFromSmarts('OCC(O)COC(=O)[C]')
    if not mol.HasSubstructMatch(acyl_glycerol_pattern):
        return False, "No glycerol backbone with acyl group at position 1 found"

    # Free hydroxyl group at position 2
    free_oh_pattern = Chem.MolFromSmarts('OC[C](O)CO[P](=O)(O)OCCN')
    if not mol.HasSubstructMatch(free_oh_pattern):
        return False, "No free hydroxyl group at position 2 found"

    return True, "Molecule is a 1-O-acylglycerophosphoethanolamine"