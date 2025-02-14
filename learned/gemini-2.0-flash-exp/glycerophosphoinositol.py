"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol has a glycerol backbone with two fatty acid chains and a phosphate group,
    where the phosphate is esterified with inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for phosphate group attached to glycerol
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4][OX2][P](=[OX1])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
          return False, "No phosphate group directly attached to glycerol"

    # 3. Check for inositol ring attached to phosphate
    inositol_pattern = Chem.MolFromSmarts("[C]1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    
    phosphorylated_inositol_pattern = Chem.MolFromSmarts("[OP](=[OX1])([OX2])-[C]1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    
    if not mol.HasSubstructMatch(phosphorylated_inositol_pattern):
        return False, "No phosphorylated inositol group found"


    # 4. Check for two fatty acid chains (long carbon chains attached to the glycerol via ester bonds)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) < 2:
       return False, f"Found {len(ester_matches)} ester groups, need at least 2 (for 2 fatty acids)"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # 5. Check rotatable bonds for fatty acid chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Fatty acid chains are too short"

    # Check for carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 15:
        return False, "Too few carbons for glycerophosphoinositol"
    
    if p_count != 1:
        return False, "Must have exactly 1 phosphate group"
        
    if o_count < 10:
        return False, "Too few oxygen atoms for glycerophosphoinositol"

    return True, "Contains glycerol backbone, 2 fatty acid chains, inositol, and a phosphate group attached to inositol and glycerol"