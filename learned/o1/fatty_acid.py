"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid derived from or contained
    in esterified form in an animal or vegetable fat, oil, or wax.
    Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched
    and even-numbered), which may be saturated or unsaturated.
    The term can also include all acyclic aliphatic carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group(s)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O-]')  # Deprotonated carboxylate
    carboxylic_acid_pattern_protonated = Chem.MolFromSmarts('C(=O)O')  # Protonated carboxylic acid
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_acid_matches += mol.GetSubstructMatches(carboxylic_acid_pattern_protonated)
    num_carboxylic_acids = len(carboxylic_acid_matches)
    if num_carboxylic_acids == 0:
        return False, "No carboxylic acid group found"

    # Check for aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    if aromatic_ring_count > 0:
        return False, "Contains aromatic rings, not typical for fatty acids"

    # Calculate total number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check that there are at least 4 carbons and up to 28 carbons
    if num_carbons < 4:
        return False, f"Too few carbon atoms ({num_carbons}), fatty acids have at least 4 carbons"
    if num_carbons > 28:
        return False, f"Too many carbon atoms ({num_carbons}), natural fatty acids have up to 28 carbons"

    # Allow functional groups like hydroxyl, methoxy, halogens (as seen in examples)
    # Ensure no uncommon elements are present
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains uncommon atom types ({atom.GetSymbol()}) not typical in fatty acids"

    return True, "Molecule is classified as a fatty acid"