"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
A sterol ester obtained by formal condensation of the carboxy group of any 
carboxylic acid with the 3-hydroxy group of cholesterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]12[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]([#6][#6][#6]2[#6][#6][#6]1)[#6]34")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for cholesterol-specific features
    # Including the double bond at C5-C6 and methyl groups
    cholesterol_pattern = Chem.MolFromSmarts("[#6]-[#6]1-[#6]-[#6]-[#6]2=[#6]-[#6]-[#6]-[#6]3-[#6]-[#6]-[#6]-[#6]-[#6]-3-[#6]-[#6]-2-[#6]-1")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "Missing characteristic cholesterol features"

    # Verify ester is attached at position 3 (the only oxygen should be part of the ester)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:  # One for ester C=O and one for ester C-O-
        return False, f"Found {oxygen_count} oxygens, should be exactly 2 for cholesteryl ester"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Cholesteryl esters are typically >500 Da
        return False, "Molecular weight too low for cholesteryl ester"

    # Count carbons - cholesteryl esters typically have >30 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 29:
        return False, "Too few carbons for cholesteryl ester"

    # Check for reasonable chain length of fatty acid part
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short for typical cholesteryl ester"

    return True, "Contains cholesterol core with ester-linked fatty acid at position 3"