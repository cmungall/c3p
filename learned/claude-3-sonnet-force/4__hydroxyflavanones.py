"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: CHEBI:28096 4'-hydroxyflavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4_hydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a flavanone core with a hydroxy substituent at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavanone core pattern (two fused rings with C=O in one ring)
    flavanone_pattern = Chem.MolFromSmarts("C1=CC(=O)C2=C(C=C1)C=CC=C2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core found"

    # Look for hydroxy group at position 4'
    hydroxy_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    hydroxy_match = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_match:
        return False, "No hydroxy group at position 4'"

    # Get atom indices of matched hydroxy group
    hydroxy_atom_idx = hydroxy_match[0][3]

    # Check if hydroxy group is connected to flavanone core at position 4'
    flavor_core_ring = mol.GetRingInfo().BondRings()[0]
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in flavor_core_ring and bond.GetEndAtomIdx() == hydroxy_atom_idx:
            if bond.GetBondType() == Chem.BondType.SINGLE:
                return True, "Contains flavanone core with hydroxy group at position 4'"
            else:
                return False, "Hydroxy group not connected via single bond"

    return False, "Hydroxy group not connected to flavanone core at position 4'"