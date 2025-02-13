"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:23066 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems have a characteristic beta-lactam ring fused to a 5-membered ring
    with specific substitution patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core carbapenem structure (beta-lactam fused to 5-membered ring)
    # [C]1[C](=O)N2[C][C][C]12 represents the basic skeleton
    carbapenem_core = Chem.MolFromSmarts("[C]1[C](=O)N2[C][C][C]12")
    if not mol.HasSubstructMatch(carbapenem_core):
        return False, "Missing carbapenem core structure (fused beta-lactam and 5-membered ring)"

    # Check for double bond in the 5-membered ring
    # This is characteristic of carbapenems
    double_bond_pattern = Chem.MolFromSmarts("[C]1[C](=O)N2[C]=C[C]12")
    if not mol.HasSubstructMatch(double_bond_pattern):
        # Also check alternative double bond position
        alt_double_bond = Chem.MolFromSmarts("[C]1[C](=O)N2C=C[C]12")
        if not mol.HasSubstructMatch(alt_double_bond):
            return False, "Missing characteristic double bond in 5-membered ring"

    # Check for carboxylic acid group (common in carbapenems)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        # Also check for carboxylate
        carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
        if not mol.HasSubstructMatch(carboxylate_pattern):
            return False, "Missing carboxylic acid/carboxylate group"

    # Check ring size (should have a 4-membered and a 5-membered ring)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not (4 in ring_sizes and 5 in ring_sizes):
        return False, "Incorrect ring sizes - needs 4 and 5 membered rings"

    # Check for common substituents often found in carbapenems
    # Sulfur-containing groups
    sulfur_pattern = Chem.MolFromSmarts("[S]")
    hydroxyethyl_pattern = Chem.MolFromSmarts("CC(O)")
    
    # Count number of nitrogens (should have at least the one in the beta-lactam)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "Missing required nitrogen atom"

    # Additional check for molecular weight (most carbapenems are between 250-650 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 800:
        return False, "Molecular weight outside typical range for carbapenems"

    return True, "Contains carbapenem core structure with characteristic substitution patterns"