"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:35755 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a single carboxylic acid group with a ketone group
    at the alpha position (2-oxo).

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

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if len(carboxylic_matches) == 0:
        return False, "No carboxylic acid group found"
    elif len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Find the carboxylic acid carbon atom
    acid_carbon = carboxylic_matches[0][0]
    
    # Pattern for 2-oxo monocarboxylic acid - more specific than before
    alpha_keto_pattern = Chem.MolFromSmarts("[#6]-[CX3](=O)-[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(alpha_keto_pattern)
    
    if not matches:
        return False, "No 2-oxo monocarboxylic acid pattern found"

    # Check that we found the same carboxylic acid carbon
    for match in matches:
        if match[3] == acid_carbon:
            # Get neighboring atoms of the alpha carbon (ketone carbon)
            alpha_carbon = match[1]
            alpha_neighbors = mol.GetAtomWithIdx(alpha_carbon).GetNeighbors()
            
            # Count number of ketone groups attached to alpha carbon
            ketone_count = 0
            for neighbor in alpha_neighbors:
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
                    ketone_count += 1
            
            if ketone_count > 1:
                return False, "Multiple ketone groups on alpha carbon"

            # Check that the ketone is not part of a larger conjugated system
            conjugated_pattern = Chem.MolFromSmarts("[CX3](=O)-[CX3](=O)-[CX3]=O")
            if mol.HasSubstructMatch(conjugated_pattern):
                return False, "Ketone is part of a conjugated system"

            # Check that the ketone carbon is not part of a ring containing the acid
            ring_info = mol.GetRingInfo()
            if ring_info.NumAtomRings(alpha_carbon) > 0 and ring_info.NumAtomRings(acid_carbon) > 0:
                # Check if both atoms are in the same ring
                for ring in ring_info.AtomRings():
                    if alpha_carbon in ring and acid_carbon in ring:
                        return False, "Ketone and acid are part of the same ring"

            return True, "Contains a valid 2-oxo monocarboxylic acid pattern"

    return False, "No valid 2-oxo monocarboxylic acid pattern found"