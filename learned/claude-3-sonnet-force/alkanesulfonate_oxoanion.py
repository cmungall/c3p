"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:52526 alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate (-SO3-) group attached to an alkane chain.

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
    
    # Look for sulfonate group pattern
    sulfonate_pattern = Chem.MolFromSmarts("[S+3(=O)(-[O-])(-[O-])]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    if len(sulfonate_matches) != 1:
        return False, f"Found {len(sulfonate_matches)} sulfonate groups, need exactly 1"
    
    # Check for alkane chain attached to sulfonate group
    sulfonate_atom_idx = sulfonate_matches[0][0]
    sulfonate_atom = mol.GetAtomWithIdx(sulfonate_atom_idx)
    
    # Check for at least one carbon neighbor
    has_alkane_chain = False
    for neighbor in sulfonate_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon
            has_alkane_chain = True
            break
    
    if not has_alkane_chain:
        return False, "No alkane chain attached to sulfonate group"
    
    # Check for overall negative charge
    charge = Chem.GetFormalCharge(mol)
    if charge >= 0:
        return False, "Molecule must have an overall negative charge"
    
    return True, "Molecule contains a sulfonate group attached to an alkane chain with overall negative charge"