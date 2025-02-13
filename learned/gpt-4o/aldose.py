"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure patterns
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)[CH2]")
    hemiacetal_pattern = Chem.MolFromSmarts("O[C@]1([C@@H](O)[C@H](O)[C@H](O)C1)")

    # Check for open-chain aldose
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Expected polyhydroxy structure (-OH groups)
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) < 2:
            return False, f"Not enough hydroxyl groups found, got {len(hydroxyl_matches)}"

        carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
        if carbon_count < 3:
            return False, f"Insufficient carbon backbone, found {carbon_count}"

        return True, "Structure is consistent with open-chain form of an aldose"

    # Check for cyclic (hemiacetal) form
    if mol.HasSubstructMatch(hemiacetal_pattern):
        return True, "Structure is consistent with cyclic hemiacetal form of an aldose"

    # In the loop above, even if other criteria appeared valid, failure to match either pattern returns False.
    return False, "Structure does not fit criteria for an aldose"