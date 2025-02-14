"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:38199 beta-lactam antibiotic
An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("[C&R]1(N[C&R]([C&R]1=[O])=[O])")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for nitrogen-containing heterocycle
    heterocycle_pattern = Chem.MolFromSmarts("*1~*~*~*~*~*1")
    heterocycle_pattern = AllChem.AddAtomProp(heterocycle_pattern, "N")
    heterocycle_matches = mol.GetSubstructMatches(heterocycle_pattern)
    if not heterocycle_matches:
        return False, "No nitrogen-containing heterocycle found"

    # Check for amide or amine groups (for antibiotic activity)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amine_pattern = Chem.MolFromSmarts("N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amide_matches and not amine_matches:
        return False, "No amide or amine groups found (required for antibiotic activity)"

    # Check for carboxylic acid group (common in beta-lactam antibiotics)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return True, "Contains beta-lactam ring, nitrogen-containing heterocycle, and amide/amine groups"
    else:
        return True, "Contains beta-lactam ring, nitrogen-containing heterocycle, amide/amine groups, and carboxylic acid group"