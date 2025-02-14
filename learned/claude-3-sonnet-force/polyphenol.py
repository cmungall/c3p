"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:36328 polyphenol

Polyphenols are defined as members of the class of phenols that contain 2 or more benzene rings,
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

    # Count benzene rings
    ring_info = mol.GetRingInfo()
    n_benzene_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring))
    if n_benzene_rings < 2:
        return False, "Less than 2 benzene rings"

    # Check for hydroxy groups on benzene rings
    hydroxy_on_benzene = False
    for ring in ring_info.AtomRings():
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
                    hydroxy_on_benzene = True
                    break
    if not hydroxy_on_benzene:
        return False, "No hydroxy groups on benzene rings"

    # Check for phenolic oxygens
    n_phenolic_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if n_phenolic_oxygens < 2:
        return False, "Less than 2 phenolic oxygens"

    return True, "Contains 2 or more benzene rings, each with at least one hydroxy group"