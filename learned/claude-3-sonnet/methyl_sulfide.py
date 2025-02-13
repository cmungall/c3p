"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:33108 methyl sulfide
Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide with at least one methyl group attached to the sulfur.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check if any sulfur atom has a methyl group attached
    has_methyl_sulfide = False
    for sulfur in sulfur_atoms:
        for neighbor in sulfur.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if carbon has 3 hydrogen neighbors (methyl group)
                if len([n for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1]) == 3:
                    has_methyl_sulfide = True
                    break
        if has_methyl_sulfide:
            break
    
    if not has_methyl_sulfide:
        return False, "No methyl group attached to sulfur"
    
    # Check if sulfide is aliphatic (no aromatic rings)
    if mol.GetAromaticRingCount() > 0:
        return False, "Contains aromatic rings, not aliphatic"
    
    return True, "Contains aliphatic sulfide with methyl group attached to sulfur"