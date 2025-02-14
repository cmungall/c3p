"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is an ester formed between a fatty acid and a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, must have exactly one"

    # Look for glycerol backbone (C-C-C with 3 oxygens attached). If present, its not a wax ester
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if mol.HasSubstructMatch(glycerol_pattern):
      return False, "Glycerol backbone detected; not a wax ester"

    # Check for fatty acid chains (long carbon chains attached to esters)
    # We need at least 3 carbons on each side of the ester
    fatty_acid_pattern1 = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[OX2]") 
    fatty_acid_pattern2 = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX3](=[OX1])")
    fatty_acid_matches1 = mol.GetSubstructMatches(fatty_acid_pattern1)
    fatty_acid_matches2 = mol.GetSubstructMatches(fatty_acid_pattern2)

    if len(fatty_acid_matches1) == 0 or len(fatty_acid_matches2) == 0 :
        return False, f"Missing long carbon chains"


    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8: # min of 8 to account for 3 carbons on each side plus one ester
      return False, "Chains too short to be fatty acid and fatty alcohol"
    
    # Check molecular weight - wax esters typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for wax ester"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
      return False, "Too few carbons for a wax ester"
    
    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
      return False, "Must have exactly 2 oxygens (1 ester group)"

    return True, "Contains one ester group linking two fatty chains"