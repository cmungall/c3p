"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:35341 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid has a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (4 fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C][C]2[C]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) connected to two carbons in the first ring
    ketone_pattern = Chem.MolFromSmarts("[CH2][C](=O)[CH2]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-beta configuration
    # In 5-beta steroids, the A/B ring junction is trans
    # We can look for the specific stereochemistry pattern
    # The pattern checks for the characteristic trans fusion at the A/B ring junction
    # with the hydrogen at position 5 being up (beta)
    ab_junction_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]2[CH2]")
    if not mol.HasSubstructMatch(ab_junction_pattern):
        return False, "Not a 5-beta configuration"

    # Additional validation: Check carbon count (steroids typically have 19-27 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"

    # If we get here, it's a 3-oxo-5beta-steroid
    return True, "Molecule contains steroid core with 3-oxo group and 5-beta configuration"