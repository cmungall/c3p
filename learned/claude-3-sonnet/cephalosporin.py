"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins contain a beta-lactam ring fused to a 6-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for the core cephalosporin structure:
    # Beta-lactam fused to 6-membered dihydrothiazine ring
    # [S] represents sulfur atom in the 6-membered ring
    # The numbers indicate the atom mapping to ensure correct connectivity
    ceph_core = Chem.MolFromSmarts("""
        [S:1]-[C:2]-[C:3]=C-[N:4]-1-C(=O)-[C:5]-[C:6]-1
        """)
    
    if not mol.HasSubstructMatch(ceph_core):
        return False, "Missing cephalosporin core structure (fused beta-lactam and dihydrothiazine rings)"

    # Check for carboxylic acid group typically present at position 2
    carboxyl = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    
    if not (mol.HasSubstructMatch(carboxyl) or mol.HasSubstructMatch(carboxylate)):
        return False, "Missing carboxylic acid/carboxylate group"

    # Check for amide group at position 7
    amide = Chem.MolFromSmarts("[NH]C(=O)")
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide group at position 7"

    # Additional check to ensure proper ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for cephalosporin structure"

    # Check for common substituents that might indicate a cephalosporin
    common_substituents = [
        (Chem.MolFromSmarts("NC(=O)"), "carboxamide"),
        (Chem.MolFromSmarts("OC(=O)"), "carboxyl"),
        (Chem.MolFromSmarts("n1nnn[nH]1"), "tetrazole"),
        (Chem.MolFromSmarts("c1csc(N)n1"), "aminothiazole"),
        (Chem.MolFromSmarts("C=C"), "vinyl"),
    ]
    
    found_substituents = []
    for pattern, name in common_substituents:
        if pattern and mol.HasSubstructMatch(pattern):
            found_substituents.append(name)

    # Final verification - molecule should have sufficient complexity
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:  # Typical cephalosporins are larger
        return False, "Molecule too small to be a cephalosporin"

    return True, f"Contains cephalosporin core structure (beta-lactam fused to dihydrothiazine ring) with typical substituents: {', '.join(found_substituents)}"