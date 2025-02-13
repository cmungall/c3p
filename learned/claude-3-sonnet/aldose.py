"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: aldose sugars
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes H[CH(OH)]nC(=O)H (n >= 2) or their cyclic hemiacetal forms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of unwanted atoms
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:  # Only H, C, O allowed
            return False, f"Contains non-sugar atoms (not H, C, or O)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an aldose"
    if o_count < 3:
        return False, "Too few oxygens for an aldose"

    # Check for aromatic rings
    if rdMolDescriptors.CalcNumAromaticRings(mol) > 0:
        return False, "Contains aromatic rings"

    # Check for carboxylic acids (exclude uronic acids)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxylic acid group"

    # Check for modified sugars
    ether_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4]")  # Non-ring ethers
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester group"
    
    # Only allow cyclic ethers (hemiacetals)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    for match in ether_matches:
        atoms = [mol.GetAtomWithIdx(i) for i in match]
        if not all(atom.IsInRing() for atom in atoms):
            return False, "Contains non-cyclic ether (modified sugar)"

    # Pattern for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    
    # Pattern for hemiacetal
    hemiacetal_pattern = Chem.MolFromSmarts("[O;R]-[CH1](-[O])-[C,O]")
    
    # Pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    
    # Count hydroxyl groups
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check ring count
    rings = mol.GetRingInfo()
    if rings.NumRings() > 1:
        return False, "Too many rings for aldose"
    
    # If cyclic, check ring size
    if rings.NumRings() == 1:
        ring_sizes = [len(r) for r in rings.AtomRings()]
        if not any(size in [5,6] for size in ring_sizes):
            return False, "Ring size not typical for aldose"

    # Check for either aldehyde or hemiacetal
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    
    if not (has_aldehyde or has_hemiacetal):
        return False, "No aldehyde or hemiacetal group found"

    # Pattern for carbon chain with hydroxyls
    polyhydroxy_chain = Chem.MolFromSmarts("[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]")
    if not mol.HasSubstructMatch(polyhydroxy_chain):
        return False, "No suitable carbon chain found"

    # Check hydroxyl ratio
    ratio_oh_to_c = hydroxyl_matches / c_count
    if ratio_oh_to_c < 0.3:
        return False, "Too few hydroxyls relative to carbons"

    if has_aldehyde:
        return True, "Contains aldehyde group with multiple hydroxyls"
    else:
        return True, "Contains cyclic hemiacetal form with multiple hydroxyls"