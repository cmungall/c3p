"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a sulfonate directly attached to an aliphatic carbon
    sulfonate_pattern = Chem.MolFromSmarts("C-S([O-])(=O)(=O)")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate group attached to an aliphatic carbon found"

    # Check for valid R group attachments
    for match in matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        
        # Ensure the carbon is not further part of an aromatic system
        if not carbon_atom.GetIsAromatic():
            # Check attached groups to the carbon
            valid_connections = any(n.GetAtomicNum() in [1, 6, 7, 8] for n in carbon_atom.GetNeighbors())
            if valid_connections:
                return True, "Found an appropriately attached sulfonate group with a diverse R group"

    return False, "Sulfonate group not properly attached to valid carbon or R group connections"