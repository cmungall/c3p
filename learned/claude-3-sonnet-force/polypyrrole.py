"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:52236 polypyrrole
A compound composed of two or more pyrrole units.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyrrole substructure
    pyrrole_pattern = Chem.MolFromSmarts("c1ccnc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Require at least two pyrrole units
    if len(pyrrole_matches) < 2:
        return False, "Less than two pyrrole units found"
    
    # Count number of nitrogen atoms
    n_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Require at least two nitrogen atoms (from pyrrole units)
    if n_nitrogens < 2:
        return False, "Less than two nitrogen atoms found"
    
    # Check for conjugated ring system
    conjugated_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        bond_orders = [mol.GetBondBetweenAtoms(ring[i], ring[i-1]).GetBondTypeAsDouble() for i in range(len(ring))]
        if sum(bond_orders) % 2 == 0 and any(order == 2 for order in bond_orders):
            conjugated_rings.append(ring)
    
    if not conjugated_rings:
        return False, "No conjugated ring system found"
    
    # Check molecular weight - polypyrroles typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for polypyrrole"
    
    return True, "Contains two or more pyrrole units in a conjugated ring system"