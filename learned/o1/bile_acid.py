"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Bile Acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid with a steroid backbone,
    hydroxyl groups at specific positions, and a carboxylic acid side chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Standardize molecule (tautomeric forms, charges, etc.)
    try:
        # Add hydrogens
        mol = Chem.AddHs(mol)
    except:
        return False, "Error processing molecule"
    
    # Define steroid backbone (four fused rings with specific ring sizes)
    steroid_pattern = Chem.MolFromSmarts("""
        [#6]1[#6][#6][#6]2[#6](=[#6]1)[#6][#6]3[#6](=[#6]2)[#6][#6]4[#6](=[#6]3)[#6][#6][#6]4
    """)
    # The above SMARTS is a simplified pattern for the steroid skeleton
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 5beta-configuration
    # In practice, stereochemistry can be complex to verify; we'll check for correct chirality centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    required_centers = {
        5: 'S'  # Carbon 5 with beta configuration (stereochemistry)
    }
    for atom_idx, chirality in chiral_centers:
        atom_number = mol.GetAtomWithIdx(atom_idx).GetAtomicNum()
        if atom_number == 6 and atom_idx in required_centers:
            if chirality != required_centers[atom_idx]:
                return False, "Incorrect 5beta-configuration"
    
    # Check for hydroxyl groups at positions 3, 7, 12 (optional)
    hydroxy_positions = [3, 7, 12]
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6:
                # Could check for specific atom indices corresponding to positions
                hydroxy_count +=1
    if hydroxy_count < 1:
        return False, "No hydroxyl groups found"
    
    # Check for carboxylic acid side chain
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid side chain found"
    
    return True, "Molecule matches bile acid structural features"