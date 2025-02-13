"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35519 beta-lactam antibiotic
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

    # Look for penam (penicillin) nucleus
    penam_pattern = Chem.MolFromSmarts("[NR]1[CR]2[CR][CR][NR]1[CR]2=O.[SX2]")
    if mol.HasSubstructMatch(penam_pattern):
        return True, "Contains a penam (penicillin) nucleus"

    # Look for cepham (cephalosporin) nucleus
    cepham_pattern = Chem.MolFromSmarts("[NR]1[CR]2[CR][CR][NR]1[CR]2=O.[SX2][CX4]")
    if mol.HasSubstructMatch(cepham_pattern):
        return True, "Contains a cepham (cephalosporin) nucleus"

    # Look for carbapenem nucleus
    carbapenem_pattern = Chem.MolFromSmarts("[NR]1[CR]2[CR][CR][NR]1[CR]2=O")
    if mol.HasSubstructMatch(carbapenem_pattern):
        return True, "Contains a carbapenem nucleus"

    # Look for monobactam nucleus
    monobactam_pattern = Chem.MolFromSmarts("[NR]1[CR]2[CR][CR][NR]1[CR]2=O.[NX3]")
    if mol.HasSubstructMatch(monobactam_pattern):
        return True, "Contains a monobactam nucleus"

    # Look for carboxyl groups and amine groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    amine_pattern = Chem.MolFromSmarts("[N;H2,H1,H0]")
    if mol.HasSubstructMatch(carboxyl_pattern) and mol.HasSubstructMatch(amine_pattern):
        return True, "Contains carboxyl and amine groups"

    return False, "No characteristic features of beta-lactam antibiotics found"