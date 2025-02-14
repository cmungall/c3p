"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a cytidine diphosphate group linked to a diacylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CDP core
    # Nucleoside with ribose and 2 phosphate
    cdp_core_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H]([C@@H](O[C@H]1[n]2[c]([n][cH][c]([n]2)[N])=O)O)COP(=O)([O])OP(=O)([O])O')
    if not mol.HasSubstructMatch(cdp_core_pattern):
       return False, "CDP core not found"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for a phosphate ester linkage to the glycerol (C-O-P)
    phosphate_ester_pattern = Chem.MolFromSmarts("[CX4][OX2][P]")
    phosphate_ester_matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    
    if len(phosphate_ester_matches) < 1:
        return False, f"Missing phosphate ester group"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - CDP-diacylglycerols typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for CDP-diacylglycerol"
    
    # Count carbons, oxygens, phosphorus and nitrogen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if c_count < 20:
        return False, "Too few carbons for CDP-diacylglycerol"
    if o_count < 14:
        return False, "Too few oxygens for CDP-diacylglycerol"
    if p_count != 2:
        return False, "Must have exactly 2 phosphorus atoms (diphosphate group)"
    if n_count != 3:
        return False, "Must have 3 nitrogen atoms in CDP base"
    
    return True, "Contains CDP core with diacylglycerol attached"