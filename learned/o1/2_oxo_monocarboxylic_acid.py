"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:37727 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is any monocarboxylic acid having a 2-oxo substituent.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for monocarboxylic acid group (protonated or deprotonated)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[$([CX3](=O)[O-]),$([CX3](=O)[OH])]')

    # Find matches for the carboxylic acid group
    acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, need exactly 1"

    # Get the index of the carboxylic acid carbon atom
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # Initialize list to store alpha carbons (carbons adjacent to the carboxylic acid group)
    alpha_carbons = []
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            alpha_carbons.append(neighbor)

    # Check if there is at least one alpha carbon
    if len(alpha_carbons) == 0:
        return False, "No alpha carbon found adjacent to carboxylic acid group"

    # Check each alpha carbon for the presence of a ketone (=O) group
    alpha_oxo_found = False
    for alpha_carbon in alpha_carbons:
        # Iterate over bonds connected to the alpha carbon
        for bond in alpha_carbon.GetBonds():
            neighbor = bond.GetOtherAtom(alpha_carbon)
            # Check for a double bond to an oxygen atom
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                alpha_oxo_found = True
                break
        if alpha_oxo_found:
            break

    if not alpha_oxo_found:
        return False, "Alpha carbon does not have a ketone (=O) group"

    return True, "Molecule is a 2-oxo monocarboxylic acid"