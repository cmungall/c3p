"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group connected to choline
    # [O,OH]-P(=O)([O-])OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("[OX2,OH][PX4](=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for glycerol backbone (C-C-C with oxygens)
    # The carbons can have different substituents (H, acyl, alkyl)
    glycerol_pattern = Chem.MolFromSmarts("[OX2,OH][CH2X4][CHX4][CH2X4][OX2,OH]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count key atoms
    n_atoms = mol.GetNumAtoms()
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)  # Phosphorus
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)   # Nitrogen
    
    if p_count != 1:
        return False, f"Must have exactly 1 phosphorus atom, found {p_count}"
    if n_count != 1:
        return False, f"Must have exactly 1 nitrogen atom, found {n_count}"

    # Check for quaternary nitrogen (choline)
    quat_n_pattern = Chem.MolFromSmarts("[NX4+]")
    if not mol.HasSubstructMatch(quat_n_pattern):
        return False, "No quaternary ammonium group found"

    # Look for common substituents
    # Ester groups
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    # Ether groups
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    
    n_esters = len(mol.GetSubstructMatches(ester_pattern))
    n_ethers = len(mol.GetSubstructMatches(ether_pattern))
    
    if n_esters == 0 and n_ethers == 0:
        return False, "No ester or ether substituents found"

    # Additional check for phosphate-glycerol connection
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[OX2,OH][PX4](=O)([O-])OC[CH]")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "Phosphate not properly connected to glycerol backbone"

    return True, "Contains glycerol backbone with phosphocholine group and appropriate substituents"