"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: cardiac glycoside
Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for steroid core (tetracyclic ring system)
    steroid_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]1)[C]1[C][C][C]3[C]([C]1[C][C]2)[C][C][C]2[C][C][C][C][C]23")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"
        
    # Check for butenolide/furanone ring
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
        
    # Check for sugar moiety (pyranose ring with multiple OH groups)
    sugar_pattern = Chem.MolFromSmarts("[OX2H1][CX4H1][CX4H1][CX4H1][CX4H1][CX4H1]O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
        
    # Check for glycosidic linkage (C-O-C between sugar and steroid)
    glycosidic_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
        
    # Count oxygen atoms (should have multiple due to OH groups and sugar)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for cardiac glycoside"
        
    # Check molecular weight (typically >500 Da for cardiac glycosides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for cardiac glycoside"
        
    # Count rings (should have multiple rings including steroid core + sugar)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:  # steroid (4) + at least 1 sugar ring
        return False, "Too few rings for cardiac glycoside"
        
    # Count hydroxyl groups (should have multiple)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 3:
        return False, "Too few hydroxyl groups"

    return True, "Contains steroid core with lactone ring and sugar residue(s)"