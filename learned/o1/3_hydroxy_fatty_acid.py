"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:25418 3-hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxy group at the beta- or 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Check each carboxylic acid group
    for match in carboxylic_acid_matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Get alpha carbon (carbon adjacent to carboxyl carbon)
        alpha_carbons = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not alpha_carbons:
            continue  # No alpha carbon found
        alpha_carbon = alpha_carbons[0]

        # Get beta carbon (carbon adjacent to alpha carbon, not carboxyl carbon)
        beta_carbons = [nbr for nbr in alpha_carbon.GetNeighbors()
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carboxyl_carbon_idx]
        if not beta_carbons:
            continue  # No beta carbon found
        beta_carbon = beta_carbons[0]

        # Check if beta carbon has a hydroxy group attached
        has_hydroxy = False
        for nbr in beta_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1:
                has_hydroxy = True
                break
        if not has_hydroxy:
            continue  # Beta carbon does not have hydroxy group

        # Check that the molecule is a fatty acid (long aliphatic chain)
        c_chain = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if num_carbons < 4:
            return False, "Carbon chain too short to be a fatty acid"

        return True, "Molecule is a 3-hydroxy fatty acid with hydroxy group at beta position"

    return False, "No hydroxy group found at beta position relative to carboxylic acid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25418',
        'name': '3-hydroxy fatty acid',
        'definition': 'Any fatty acid with a hydroxy functional group in the beta- or 3-position. '
                      'beta-Hydroxy fatty acids accumulate during cardiac hypoxia, and can also be '
                      'used as chemical markers of bacterial endotoxins.',
        'parents': ['CHEBI:35366', 'CHEBI:15846']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8
    }
}