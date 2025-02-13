"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:24851 polyphenol

A polyphenol is defined as a member of the class of phenols that contain 2 or more benzene rings
each of which is substituted by at least one hydroxy group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count aromatic rings
    aromatic_rings = mol.GetRingInfo().AtomRings()
    num_aromatic_rings = sum(1 for ring in aromatic_rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if num_aromatic_rings < 2:
        return False, "Less than 2 aromatic rings"
    
    # Check for hydroxyl groups on aromatic rings
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    has_hydroxy_on_aromatic = False
    for atom_idx in aromatic_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            has_hydroxy_on_aromatic = True
            break
    if not has_hydroxy_on_aromatic:
        return False, "No hydroxyl groups on aromatic rings"
    
    # Check if all aromatic rings have at least one hydroxyl group
    ring_with_hydroxy = set()
    for atom_idx in aromatic_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            ring_idx = mol.GetAtomRingInfo().IsAtomInRingOfSize(atom_idx, 6)
            if ring_idx >= 0:
                ring_with_hydroxy.add(ring_idx)
    if len(ring_with_hydroxy) < num_aromatic_rings:
        return False, "Not all aromatic rings have a hydroxyl group"
    
    return True, "Contains 2 or more aromatic rings, each with at least one hydroxyl group"