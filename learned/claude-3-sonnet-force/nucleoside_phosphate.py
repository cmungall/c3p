"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:35489 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecular entity that is a nucleoside
    in which one or more of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nucleobase pattern
    nucleobase_pattern = Chem.MolFromSmarts("[nr3]1[nr3][nr3][nr3][nr3]1")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"
    
    # Look for sugar ring pattern
    sugar_ring_pattern = Chem.MolFromSmarts("[OX2r3][CX4r3][CX4r3][CX4r3][CX4r3][OX2r3]")
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar ring found"
    
    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"
    
    # Check if phosphate is attached to sugar hydroxy
    sugar_oh_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 and atom.IsInRing()]
    phosphate_attached = False
    for phosphate_match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(phosphate_match[0])
        for neighbor in phosphate_atom.GetNeighbors():
            if neighbor.GetIdx() in [atom.GetIdx() for atom in sugar_oh_atoms]:
                phosphate_attached = True
                break
    if not phosphate_attached:
        return False, "Phosphate not attached to sugar hydroxy"
    
    return True, "Contains a nucleobase, sugar ring, and one or more phosphate groups attached to sugar hydroxy"