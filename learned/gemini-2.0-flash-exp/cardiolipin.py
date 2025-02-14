"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two molecules of phosphatidic acid
    covalently linked to a molecule of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Central glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No central glycerol backbone found"
    
    # Search for 2 phosphate groups (each connected to 2 oxygens, one of which is part of glycerol)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])") #X2 to ensure connected to two other groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 2"

    # Search for 4 ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 4"

    # Check for 4 fatty acid chains (long carbon chains) attached to esters - at least 3 C per chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # at least 3 carbons chain
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 4:
         return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
        
    # Check for total number of glycerol carbons: should be at least 3 per phosphate group, plus 3 from the central glycerol = 9.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 9: # at least 3 glycerol carbon each + the one from central glycerol
        return False, f"Too few carbons: {c_count}, less than expected"

    # Check number of phosphorous atoms (should be 2)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 2:
        return False, f"Must have exactly 2 phosphorus atoms, found {p_count}"
    
    # Check for the number of oxygens (should be more than 12)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 12:
        return False, f"Too few oxygens: {o_count}, less than expected"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - cardiolipins typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for cardiolipin"

    return True, "Contains two phosphatidic acids linked to a glycerol backbone"