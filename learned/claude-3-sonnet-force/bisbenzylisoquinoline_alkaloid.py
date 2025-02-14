"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    These compounds contain two benzylisoquinoline units linked by ether bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic molecular properties
    # These compounds typically have:
    # - Multiple rings (>6)
    # - At least 2 nitrogen atoms
    # - Multiple oxygen atoms (for ether bridges and substituents)
    rings = rdMolDescriptors.CalcNumRings(mol)
    if rings < 6:
        return False, f"Too few rings ({rings}), need at least 6"

    # Count N and O atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2:
        return False, f"Need at least 2 nitrogen atoms, found {n_count}"
    if o_count < 2:
        return False, f"Need at least 2 oxygen atoms for ether bridges, found {o_count}"

    # Look for isoquinoline core structure (can be partially saturated)
    isoquinoline_pattern = Chem.MolFromSmarts("c1cccc2c1CCN=C2")
    tetrahydroisoquinoline_pattern = Chem.MolFromSmarts("C1CCNCc2ccccc12")
    
    iso_matches = len(mol.GetSubstructMatches(isoquinoline_pattern))
    tetra_matches = len(mol.GetSubstructMatches(tetrahydroisoquinoline_pattern))
    total_iso_units = iso_matches + tetra_matches
    
    if total_iso_units < 2:
        return False, f"Need at least 2 isoquinoline units, found {total_iso_units}"

    # Look for ether bridges (O connecting two aromatic rings)
    ether_bridge_pattern = Chem.MolFromSmarts("c-O-c")
    ether_matches = len(mol.GetSubstructMatches(ether_bridge_pattern))
    
    if ether_matches < 1:
        return False, "No ether bridges found between aromatic rings"

    # Look for common substituents
    methoxy_pattern = Chem.MolFromSmarts("cOC")
    hydroxy_pattern = Chem.MolFromSmarts("cO[H]")
    methylenedioxy_pattern = Chem.MolFromSmarts("OCO")
    
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))
    methylenedioxy_count = len(mol.GetSubstructMatches(methylenedioxy_pattern))
    
    total_substituents = methoxy_count + hydroxy_count + methylenedioxy_count
    if total_substituents < 1:
        return False, "Missing typical oxygen-containing substituents"

    # Check for N-methylation pattern (common but not required)
    n_methyl_pattern = Chem.MolFromSmarts("CN")
    n_methyl_count = len(mol.GetSubstructMatches(n_methyl_pattern))

    # Molecular weight check (these are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for bisbenzylisoquinoline alkaloid"

    reason = (f"Contains {total_iso_units} isoquinoline units, {ether_matches} ether bridges, "
             f"and appropriate substituents (MW: {mol_wt:.1f})")
    
    return True, reason