"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll molecules
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins with a fifth ring and various side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom (including different charge states)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if not mg_atoms:
        return False, "No magnesium atom found"
    if len(mg_atoms) > 1:
        return False, "More than one magnesium atom found"

    # Check for 4 nitrogens in the molecule
    n_atoms = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'])
    if n_atoms != 4:
        return False, f"Must have exactly 4 nitrogens, found {n_atoms}"

    # Simplified porphyrin-like core pattern
    porphyrin_pattern = Chem.MolFromSmarts("[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin-like core found"

    # Check for fifth ring (cyclopentanone or similar)
    fifth_ring_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#6]1")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No characteristic fifth ring found"

    # Count rings
    rings = mol.GetRingInfo()
    if rings.NumRings() < 5:
        return False, "Too few rings for chlorophyll structure"

    # Check for characteristic substituents
    substituents = []
    
    # Vinyl group (-CH=CH2)
    vinyl_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(vinyl_pattern):
        substituents.append("vinyl")
        
    # Ethyl group
    ethyl_pattern = Chem.MolFromSmarts("CC")
    if mol.HasSubstructMatch(ethyl_pattern):
        substituents.append("ethyl")
    
    # Carboxyl/ester groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(carboxyl_pattern):
        substituents.append("carboxyl/ester")

    # Long chain (phytol or similar)
    long_chain = Chem.MolFromSmarts("CCCC")
    if mol.HasSubstructMatch(long_chain):
        substituents.append("alkyl chain")

    if not substituents:
        return False, "Missing characteristic substituents"

    # Count carbons and check molecular complexity
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for chlorophyll structure"

    # Additional check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated system"

    return True, f"Contains magnesium-coordinated porphyrin core with fifth ring and {', '.join(substituents)}"