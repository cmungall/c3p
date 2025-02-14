"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:35881 secondary ammonium ion
An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for secondary ammonium ion pattern
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NH2+][CX4,CX3]([#6])[CX4,CX3]([#6])[#6]")
    matches = mol.GetSubstructMatches(secondary_ammonium_pattern)

    if not matches:
        return False, "No secondary ammonium ion pattern found"

    # Check for other heteroatoms connected to nitrogen
    for match in matches:
        try:
            n_atom = mol.GetAtomWithIdx(match[0])
            if sum(1 for atom in n_atom.GetNeighbors() if atom.GetAtomicNum() not in [1, 6, 7]) > 0:
                return False, "Secondary ammonium ion nitrogen connected to non-C/N/O atoms"
        except IndexError:
            return False, "Invalid atom index encountered"

    # Check if nitrogen is part of a ring
    is_in_ring = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in matches[0])

    # Check formal charge of nitrogen
    n_charge = sum(mol.GetAtomWithIdx(idx).GetFormalCharge() for idx in matches[0])
    if n_charge != 1:
        return False, "Nitrogen does not have a formal charge of +1"

    # Classify as secondary ammonium ion
    reason = "Contains a secondary ammonium ion group"
    if is_in_ring:
        reason += " as part of a ring"
    return True, reason