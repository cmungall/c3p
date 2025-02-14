"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid that contains three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for terminal carboxylic acid group pattern
        carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
        if not mol.HasSubstructMatch(carboxylic_acid_pattern):
            return False, "No terminal carboxylic acid group found"

        # Identify all C=C double bonds
        double_bond_pattern = Chem.MolFromSmarts("C=C")
        double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
        
        # Count non-conjugated double bonds to avoid counting within rings and branches
        unique_double_bonds = set()
        for match in double_bond_matches:
            if mol.GetAtomWithIdx(match[0]).GetDegree() < 3 and mol.GetAtomWithIdx(match[1]).GetDegree() < 3:
                # Only count double bonds part of the main chain
                unique_double_bonds.add(match)

        if len(unique_double_bonds) != 3:
            return False, f"Found {len(unique_double_bonds)} distinct double bonds along main chain, need exactly 3"

        # Check for acyclic structure
        if mol.GetRingInfo().NumRings() > 0:
            return False, "Fatty acids should be acyclic"

        # Ensure sufficient length for a typical fatty acid
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count < 12:
            return False, "Carbon chain too short for a typical fatty acid"

        return True, "Contains an acyclic fatty acid structure with exactly three distinct double bonds"
    
    except Exception as e:
        return None, f"Error occurred: {str(e)}"