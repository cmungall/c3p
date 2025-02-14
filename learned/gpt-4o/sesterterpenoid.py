"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid typically has a C25 skeleton derived from a sesterterpene,
    and may include rearrangements or modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the carbon count is around 25, which is typical for sesterterpenoid
    if c_count < 20 or c_count > 30:
        return False, "Carbon count not typical for sesterterpenoids"
    
    # Identify isoprene units pattern in terpenoids
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")  # Isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Insufficient isoprene unit matches found"
    
    # Look for common functional group patterns in terpenoids
    # For example, check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check for at least one hydroxyl group (common in terpenoids)
    if len(hydroxyl_matches) == 0:
        return False, "No typical terpenoid functional groups like hydroxyl found"
    
    # Additional validation could include patterns for rearranged carbon structures
    # Since molecular alterations and characteristic motifs are possible in sesterterpenoids

    return True, "Molecule matches basic structure with isoprene units and functional group typical for sesterterpenoids"