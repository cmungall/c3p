"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxy group on the alpha carbon
    to one of the carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Find all carboxylic acid groups
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_matches) < 2:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups, need at least 2"

    found_alpha_hydroxy = False

    # Check each carboxylic acid group
    for group in carboxylic_matches:
        carboxylic_carbon_idx = group[0]
        
        # Look for alpha carbon with hydroxy
        carboxylic_carbon = mol.GetAtomWithIdx(carboxylic_carbon_idx)
        for neighbor in carboxylic_carbon.GetNeighbors():
            if neighbor.GetSymbol() == 'C':  # Check for carbon
                # Look for hydroxy on this alpha carbon
                alpha_carbon_idx = neighbor.GetIdx()
                for alpha_neighbor in neighbor.GetNeighbors():
                    if alpha_neighbor.GetSymbol() == 'O' and alpha_neighbor.GetTotalValence() == 2:  # Check for single-bonded oxygen (hydroxyl)
                        found_alpha_hydroxy = True
                        break
            if found_alpha_hydroxy:
                break
        if found_alpha_hydroxy:
            break

    if found_alpha_hydroxy:
        return True, "Contains two carboxylic groups with a hydroxy group on an alpha carbon"
    else:
        return False, "Does not have a hydroxy group on an alpha carbon to a carboxylic acid group"