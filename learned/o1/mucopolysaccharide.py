"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Adjusted to detect molecules containing sulfur-containing heterocycles (e.g., thiazole, thiazoline)
    and disulfide bridges, as these features are present in the provided examples.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as mucopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for sulfur-containing heterocycles
    thiazole_pattern = Chem.MolFromSmarts('c1scn[cH]1')  # Thiazole ring
    thiazoline_pattern = Chem.MolFromSmarts('c1sccn1')   # Thiazoline ring
    oxazole_pattern = Chem.MolFromSmarts('c1ocn[cH]1')   # Oxazole ring
    oxazoline_pattern = Chem.MolFromSmarts('c1occn1')    # Oxazoline ring

    # Define pattern for disulfide bridge
    disulfide_pattern = Chem.MolFromSmarts('S-S')        # Disulfide bond

    # Check for sulfur-containing heterocycles
    has_thiazole = mol.HasSubstructMatch(thiazole_pattern)
    has_thiazoline = mol.HasSubstructMatch(thiazoline_pattern)
    has_oxazole = mol.HasSubstructMatch(oxazole_pattern)
    has_oxazoline = mol.HasSubstructMatch(oxazoline_pattern)
    has_disulfide = mol.HasSubstructMatch(disulfide_pattern)

    if any([has_thiazole, has_thiazoline, has_oxazole, has_oxazoline, has_disulfide]):
        reasons = []
        if has_thiazole:
            reasons.append("Contains a thiazole ring")
        if has_thiazoline:
            reasons.append("Contains a thiazoline ring")
        if has_oxazole:
            reasons.append("Contains an oxazole ring")
        if has_oxazoline:
            reasons.append("Contains an oxazoline ring")
        if has_disulfide:
            reasons.append("Contains a disulfide bridge")
        reason_str = "; ".join(reasons)
        return True, reason_str

    return False, "Does not contain sulfur-containing heterocycles or disulfide bridges"