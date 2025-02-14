"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:18115 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is a glycerophosphoglycerol, where glycerol has one of its primary
    hydroxy groups replaced by a phosphatidyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol-phosphate moiety ([CH2X4][CHX4][CH2X4]OP(=O)(O)O)
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol-phosphate moiety found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found fewer than 2 ester groups"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"
    
    # Check molecular weight - phosphatidylglycerols typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return None, "Molecular weight may be too low for phosphatidylglycerol"
    
    # Count carbons, oxygens, and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidylglycerol"
    if o_count < 8:
        return False, "Too few oxygens for phosphatidylglycerol"
    if p_count < 1:
        return False, "No phosphorus atom found"
    
    # Rule out phosphatidylcholines (PCs) and phosphatidylethanolamines (PEs)
    pc_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)")
    pe_pattern = Chem.MolFromSmarts("OCCN")
    if mol.HasSubstructMatch(pc_pattern) or mol.HasSubstructMatch(pe_pattern):
        return False, "Molecule appears to be a phosphatidylcholine or phosphatidylethanolamine"
    
    return True, "Contains glycerol-phosphate moiety with 2 fatty acid chains"