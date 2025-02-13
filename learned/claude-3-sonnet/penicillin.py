"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins have a characteristic penam core (beta-lactam fused to thiazolidine ring)
    with specific substituents: two methyls at position 2, carboxylate at position 3,
    and carboxamido at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for penam core structure (beta-lactam fused to thiazolidine ring)
    # [S]-[CH]-[N]-C(=O)-[CH]-[CH]-[N] forms the core scaffold
    penam_core = Chem.MolFromSmarts('[S]1[CH]2[N]C(=O)[CH][CH]2[N]1')
    if not mol.HasSubstructMatch(penam_core):
        return False, "No penam core structure found"

    # Check for two methyl groups at position 2
    two_methyls = Chem.MolFromSmarts('[CH3][C]([CH3])(S1)[C]')
    if not mol.HasSubstructMatch(two_methyls):
        return False, "Missing two methyl groups at position 2"

    # Check for carboxylate group at position 3
    carboxylate = Chem.MolFromSmarts('[C](S1)([CH]N)C(=O)[OH,O-]')
    if not mol.HasSubstructMatch(carboxylate):
        return False, "Missing carboxylate group at position 3"

    # Check for carboxamido group at position 6
    carboxamido = Chem.MolFromSmarts('[CH](NC(=O)*)C(=O)N')
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Missing carboxamido group at position 6"

    # Additional check for correct stereochemistry
    # We look for the characteristic trans junction between the beta-lactam and thiazolidine rings
    stereochem_pattern = Chem.MolFromSmarts('[S][C@H]1[N]C(=O)[C@@H]')
    if not mol.HasSubstructMatch(stereochem_pattern):
        return False, "Incorrect stereochemistry at ring junction"

    # Verify basic penicillin formula by counting key atoms
    # Should have exactly one sulfur atom
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count != 1:
        return False, f"Incorrect number of sulfur atoms (found {s_count}, expected 1)"

    # Should have at least 2 nitrogen atoms (one in beta-lactam, one in carboxamido)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Too few nitrogen atoms (found {n_count}, expected at least 2)"

    return True, "Contains penam core with correct substituents and stereochemistry"