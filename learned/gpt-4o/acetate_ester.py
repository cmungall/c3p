"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester, specifically a carboxylic ester where
    the carboxylic acid component is acetic acid, based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # SMARTS pattern to find the acetate group
    # Acetate group (CH3COO-): carbonyl C(=O) bonded to O and methyl group [CH3]
    # This pattern assumes [-C(=O)OCH3] where the methyl is directly attached to the ester oxygen
    acetate_pattern = Chem.MolFromSmarts("[C](=O)O[C]") 

    # Look for acetate group
    if not mol.HasSubstructMatch(acetate_pattern):
        return (False, "No acetate ester group found")
    
    # Check each match to ensure the pattern fully represents an acetate ester part of acetic acid
    matches = mol.GetSubstructMatches(acetate_pattern)
    for match in matches:
        # Verify the carbon bonded to the ester oxygen is part of a methyl group
        ester_oxygen_idx = match[2]
        connected_atoms = [neighb for neighb in mol.GetAtomWithIdx(ester_oxygen_idx).GetNeighbors()]
        if any((atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3) for atom in connected_atoms):
            return (True, "Contains acetate ester group as part of acetic acid")

    return (False, "Acetate ester group not part of acetic acid component")