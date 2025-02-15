"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:25859 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene (C30 skeleton) and may have modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polycyclic structure (at least 3 rings, as some triterpenoids may have fewer rings due to modifications)
    n_rings = mol.GetRingInfo().NumRings()
    if n_rings < 3:
        return False, f"Too few rings ({n_rings}) for a triterpenoid"

    # Check for common functional groups in triterpenoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_hydroxyl or has_carbonyl or has_carboxyl):
        return False, "No common functional groups (hydroxyl, carbonyl, carboxyl) found"

    # Check molecular weight (triterpenoids typically >300 Da, allowing for smaller derivatives)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a triterpenoid"

    # Check for long carbon chains or complex branching
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too few rotatable bonds for a triterpenoid"

    # Check for a triterpenoid-like skeleton (polycyclic with at least 20 carbons, allowing for modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a triterpenoid-like skeleton"

    return True, "Contains polycyclic structure with functional groups typical of triterpenoids"