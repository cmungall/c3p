"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerol backbone with a phosphate group attached to the
    third carbon, a choline head group and two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Glycerol Backbone with chirality
    glycerol_pattern = Chem.MolFromSmarts("[C@H]([OX2])([OX2])[CH2X4][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct chirality found"
    
    # 2. Phosphate Group attached to glycerol
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4][OX2][PX4](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to the glycerol backbone"
    
    # 3. Choline Head Group connected to phosphate
    choline_pattern = Chem.MolFromSmarts("[PX4]([OX2])([OX1])([OX2])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
         return False, "No choline head group found"
    
    #4. Two ester groups attached to the glycerol at positions 1 and 2
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2: #there are at least 2 ester bonds, but they might not be in the correct location
         return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    # ensure they are attached to the glycerol
    ester_glycerol_pattern = Chem.MolFromSmarts("[CH]([OX2][CX3](=[OX1]))([OX2][CX3](=[OX1]))[CH2X4][OX2][PX4](=[OX1])[OX2]")
    ester_glycerol_matches = mol.GetSubstructMatches(ester_glycerol_pattern)
    if len(ester_glycerol_matches) !=1:
        return False, "Ester groups are not on the 1 and 2 positions of the glycerol backbone"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
         return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6: #2 for the chain, 2 for the esters, 2 for phosphate and choline
        return False, "Chains too short to be fatty acids"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 10:
         return False, "Too few carbons for a phosphatidylcholine"

    if o_count < 7:
         return False, "Must have at least 7 oxygens"
    
    if n_count != 1:
        return False, "Must have exactly one nitrogen (choline)"
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus (phosphate)"

    return True, "Meets all criteria for phosphatidylcholine"