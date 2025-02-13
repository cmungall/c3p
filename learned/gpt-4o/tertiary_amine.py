"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three hydrocarbyl groups, 
    i.e., three carbon atoms, directly or through a simple chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over each atom in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetDegree() == 3:
            # Check if the nitrogen is bonded to 3 carbon atoms
            carbon_bonds = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    carbon_bonds += 1
                # Ensure they are fully hydrocarbyl groups (not part of amide, etc.)
                if neighbor.IsInRing() or neighbor.GetIsAromatic():
                    carbon_bonds = -1  # Invalid setting for hydrocarbyl group
                    break

            if carbon_bonds == 3:
                return True, "Contains a tertiary amine group (N bonded to 3 carbons or equivalent hydrocarbyl groups)"

    return False, "No tertiary amine group found"

# Example usages for validation
print(is_tertiary_amine("CCN(CC)CC"))  # Expected: True, "Contains a tertiary amine group (N bonded to 3 carbons)"
print(is_tertiary_amine("CCO"))       # Expected: False, "No tertiary amine group found"
print(is_tertiary_amine("CN(C)C"))    # Expected: True, "Contains a tertiary amine group (N bonded to 3 carbons)"