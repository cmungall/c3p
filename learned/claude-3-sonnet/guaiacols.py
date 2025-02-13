"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:27090 guaiacols
Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacol(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenol ring (aromatic ring with -OH substituent)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_matches:
        return False, "No phenol moiety found"
    
    # Look for methoxy (-O-CH3) substituent at ortho position
    methoxy_pattern = Chem.MolFromSmarts("Oc1ccccc1OC")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if not methoxy_matches:
        return False, "No methoxy group at ortho position to phenol ring"
    
    # Check if there is only one methoxy and one hydroxyl group
    atoms = mol.GetAtoms()
    o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (1 methoxy, 1 hydroxyl)"
    
    return True, "Contains a phenol ring with a methoxy substituent at the ortho position"