"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27279 vitamin D
Any member of a group of fat-soluble hydroxy seco-steroids that exhibit biological activity against vitamin D deficiency.
Vitamin D can be obtained from sun exposure, food and supplements and is biologically inactive and converted into
the biologically active calcitriol via double hydroxylation in the body.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D is a group of fat-soluble hydroxy seco-steroids with a specific structural pattern.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for vitamin D ring system pattern
    vit_d_pattern = Chem.MolFromSmarts("[CX4]1(CC[CX3]2[CX3]1(CCCC2=C)[H])[H]")
    if not mol.HasSubstructMatch(vit_d_pattern):
        return False, "No vitamin D ring system found"
    
    # Look for cis double bond pattern
    cis_db_pattern = Chem.MolFromSmarts("/C=C/C=C/")
    if not mol.HasSubstructMatch(cis_db_pattern):
        return False, "Missing cis double bond pattern"
    
    # Look for hydroxyl group pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found, vitamin D must have at least one"
    
    # Check for typical molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
        return False, "Molecular weight outside typical range for vitamin D"
    
    # Count carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25 or o_count < 2:
        return False, "Insufficient carbon or oxygen atoms for vitamin D"
    
    return True, "Contains vitamin D ring system, cis double bonds, and hydroxyl groups"