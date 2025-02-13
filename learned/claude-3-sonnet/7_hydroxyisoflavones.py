"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
A hydroxyisoflavone compound having a hydroxy group at the 7-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Generate 2D coordinates for the molecule
    AllChem.Compute2DCoords(mol)
    
    # Basic isoflavone core pattern:
    # [#6]1=[#6]-[#6]2=[#6](-[#6](=[#8])-[#6](=[#6]1)-[#8]-[#6]2=[#6])-[#6]3=[#6]=[#6]=[#6]=[#6]=[#6]3
    isoflavone_core = Chem.MolFromSmarts("[#6]1=[#6]-[#6]2=[#6](-[#6](=[#8])-[#6](=[#6]1)-[#8]-[#6]2=[#6])-[#6]3=[#6]=[#6]=[#6]=[#6]=[#6]3")
    
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "No isoflavone core structure found"
    
    # Pattern for 7-hydroxy position in isoflavone
    # Note: The numbering in isoflavones starts from the oxygen atom in the chromone ring
    # Position 7 is on the benzene ring of the chromone system
    hydroxy_7_pattern = Chem.MolFromSmarts("O-[c]1[c][c]2[o][c][c](-[c]3[c][c][c][c][c]3)[c](=O)[c]2[c][c]1")
    
    if not mol.HasSubstructMatch(hydroxy_7_pattern):
        return False, "No hydroxy group at position 7"
    
    # Additional checks to confirm it's specifically a 7-hydroxyisoflavone
    
    # Count the basic elements that should be present
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 3:  # Minimum 3 oxygens (chromone oxygen + carbonyl + 7-OH)
        return False, "Insufficient number of oxygen atoms"
    
    # Verify aromatic nature
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if num_aromatic_atoms < 12:  # Minimum for isoflavone core
        return False, "Insufficient aromatic system"
    
    # Check for the carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"
    
    return True, "Molecule contains isoflavone core with a hydroxy group at position 7"