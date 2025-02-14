"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is a tetrahydrofuran ring with a carbonyl group directly attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for THF ring pattern (C-C-C-C-O)
    thf_pattern = Chem.MolFromSmarts("C1CCOC1")
    if not mol.HasSubstructMatch(thf_pattern):
        return False, "No tetrahydrofuran ring found"
    
    # Check for carbonyl group directly attached to the ring (R-C(=O)-)
    carbonyl_pattern = Chem.MolFromSmarts("[C;R1](=[O])")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    if not carbonyl_matches:
          return False, "No carbonyl group attached to ring"

    # Check if the carbonyl is attached to THF ring
    thf_match = mol.GetSubstructMatch(thf_pattern)
    
    carbonyl_ring_atoms = []
    for match in carbonyl_matches:
        # Get the ring atom attached to the carbonyl
        for atom in mol.GetAtomWithIdx(match[0]).GetNeighbors():
            if atom.IsInRing():
                carbonyl_ring_atoms.append(atom.GetIdx())

    is_attached_to_thf = False

    for carbonyl_atom in carbonyl_ring_atoms:
         if carbonyl_atom in thf_match:
            is_attached_to_thf = True
            break

    if not is_attached_to_thf:
        return False, "Carbonyl group not directly attached to the tetrahydrofuran ring"
            
    return True, "Contains a tetrahydrofuran ring with a carbonyl group directly attached"