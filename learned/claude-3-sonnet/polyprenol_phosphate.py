"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal allylic hydroxy group 
of a polyprenol with 1 mol eq. of phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    diphosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_diphosphate = mol.HasSubstructMatch(diphosphate_pattern)
    
    if not (has_phosphate or has_diphosphate):
        return False, "No phosphate/diphosphate group found"

    # Look for isoprene units with flexible matching for cis/trans configurations
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-[C]=[C]-[CH2]-[CH2]")
    isoprene_pattern2 = Chem.MolFromSmarts("[CH3]-[C](-[CH3])=[C]-[CH2]")
    
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    isoprene_matches2 = len(mol.GetSubstructMatches(isoprene_pattern2))
    total_isoprene = max(isoprene_matches, isoprene_matches2)
    
    if total_isoprene < 1:
        return False, "No isoprene units found"
    
    # Check for allylic oxygen-phosphorus connection with flexible pattern
    allylic_phosphate = Chem.MolFromSmarts("[CH2,CH3]-[C]=[C]-[CH2]-[OX2]P")
    if not mol.HasSubstructMatch(allylic_phosphate):
        return False, "No phosphate group connected to allylic position"
        
    # Count carbons to verify chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:  # Minimum 1 isoprene unit (C5)
        return False, "Carbon chain too short for polyprenol"
        
    # Look for characteristic branching pattern of polyprenols
    methyl_branches = Chem.MolFromSmarts("[CH3]-[C]=[C,CH1]")
    branch_count = len(mol.GetSubstructMatches(methyl_branches))
    
    if branch_count < 1:
        return False, "No methyl branches found"

    # Check for phosphate/diphosphate at terminal position
    terminal_phosphate = Chem.MolFromSmarts("[CH2]-[OX2]P(=O)([OX2])[OX2]")
    terminal_diphosphate = Chem.MolFromSmarts("[CH2]-[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    
    if not (mol.HasSubstructMatch(terminal_phosphate) or mol.HasSubstructMatch(terminal_diphosphate)):
        return False, "Phosphate group not at terminal position"

    # Success case - provide detailed reason
    phosphate_type = "diphosphate" if has_diphosphate else "phosphate"
    return True, f"Contains polyprenol chain with {phosphate_type} group at terminal allylic position"