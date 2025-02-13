"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic molecular properties
    # Cucurbitacins are large molecules with specific atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:  # Cucurbitacins typically have 30+ carbons
        return False, "Too few carbons for cucurbitacin"
    
    if o_count < 2:  # Cucurbitacins have multiple oxygen-containing groups
        return False, "Too few oxygens for cucurbitacin"

    # Check for tetracyclic core structure (four connected rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Missing tetracyclic core structure"

    # Look for characteristic cucurbitacin substructures
    # Core tetracyclic structure with specific substitution pattern
    tetracyclic_core = Chem.MolFromSmarts("C1CC2CCC3C4CCC(C4)C3C2C1")
    if not mol.HasSubstructMatch(tetracyclic_core):
        return False, "Missing characteristic tetracyclic core"

    # Check for typical functional groups found in cucurbitacins
    
    # Ketone groups
    ketone_pattern = Chem.MolFromSmarts("[CH2,CH1,CH0]-C(=O)-[CH2,CH1,CH0]")
    ketone_matches = len(mol.GetSubstructMatches(ketone_pattern))
    
    # Hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CH2,CH1,CH0]-[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if ketone_matches == 0 and hydroxyl_matches == 0:
        return False, "Missing characteristic oxygen-containing functional groups"

    # Calculate molecular weight - cucurbitacins are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Cucurbitacins typically >450 Da
        return False, "Molecular weight too low for cucurbitacin"

    # Count rings of size 6 (characteristic of cucurbitacins)
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    six_membered_rings = sum(1 for size in ring_sizes if size == 6)
    
    if six_membered_rings < 3:
        return False, "Insufficient number of 6-membered rings"

    # Look for characteristic side chain
    side_chain = Chem.MolFromSmarts("CC(C)(O)CC")
    if not mol.HasSubstructMatch(side_chain):
        side_chain2 = Chem.MolFromSmarts("CC(C)(O)C=C")
        if not mol.HasSubstructMatch(side_chain2):
            return False, "Missing characteristic side chain"

    # If all checks pass, it's likely a cucurbitacin
    return True, "Contains tetracyclic core and characteristic cucurbitacin features"