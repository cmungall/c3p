"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:35681 secondary alcohol
A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom
which has two other carbon atoms attached to it.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all alcohol groups (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    # Check each alcohol group
    for match in alcohol_matches:
        # Get the carbon atom the -OH is attached to
        carbon_idx = mol.GetAtomWithIdx(list(mol.GetAtomWithIdx(match).GetNeighbors())[0].GetIdx())
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Check if carbon atom is saturated (no double bonds) and has 2 other carbon neighbors
        if carbon_atom.GetHybridization() == Chem.HybridizationType.SP3 and sum(1 for n in carbon_atom.GetNeighbors() if n.GetAtomicNum() == 6) == 2:
            return True, "Molecule contains a secondary alcohol group"

    return False, "No secondary alcohol group found"