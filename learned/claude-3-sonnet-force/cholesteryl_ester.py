"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: CHEBI:38048 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a sterol ester obtained by formal condensation of the carboxy group
    of any carboxylic acid with the 3-hydroxy group of cholesterol.

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

    # Look for cholesterol backbone
    cholesterol_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@@]2([C@@H]3CC[C@@]4(C)[C@H](CC[C@@H]4[C@H]3CC=C2C1)[C@@H](C)CCCC(C)C)C")
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "No cholesterol backbone found"
        
    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for acyl chain (carbon chain attached to ester)
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl chain found"

    # Count rotatable bonds to verify long acyl chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short"

    # Check molecular weight - cholesteryl esters typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for cholesteryl ester"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for cholesteryl ester"
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (ester group)"

    return True, "Contains cholesterol backbone with acyl chain attached via ester bond"