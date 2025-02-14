"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins have a flavan-3-ol skeleton and may have various hydroxyl and other substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Catechins must have a flavan-3-ol skeleton, which has a specific bicyclic core
    # SMARTS pattern for the flavan-3-ol backbone
    # [C@H]1, [C@@H]1 denote chiral centers, which may be relevant
    flavan3ol_pattern = Chem.MolFromSmarts("C1OC2=C(C=CC=C2O)CC1O")  # This is a generic pattern and may not match all specific catechin structures
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol skeleton found"

    # Check for hydroxyl groups which are typical in catechins
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if hydroxyl_count < 2:  # Typically more, but let's set at least 2 for basic confirmation
        return False, f"Too few hydroxyl groups for a catechin, found {hydroxyl_count}"

    return True, "Molecule matches flavan-3-ol skeleton and has sufficient hydroxyl groups"