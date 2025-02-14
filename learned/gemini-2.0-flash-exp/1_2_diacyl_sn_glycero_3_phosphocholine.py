"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone with correct stereochemistry
    #   The 'sn' stereochemistry is enforced by specifying the chiral center
    glycerol_pattern = Chem.MolFromSmarts("[C@H]([OX2])[CH2X4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
         return False, "No glycerol backbone with correct stereochemistry found"
        
    # 2. Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
       return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # 3. Look for the phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX1-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
         return False, "No phosphocholine group found"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
         return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    #4. Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"
    
    # 5. Count key atoms P, N, O.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if p_count != 1:
       return False, "Must have exactly 1 phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen atom"
    if o_count < 7: # There are 2 in glycerol, 2 in each ester, and 3 in phosphate
        return False, f"Must have at least 7 oxygens, found {o_count}"
    
    
    # 6. Check for negative charge on the phosphate group
    has_neg_charge = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
           has_neg_charge = True
           break

    if not has_neg_charge:
       return False, "Phosphate group must have a negative charge"
    

    return True, "Contains a glycerol backbone, two fatty acid chains attached via ester bonds, and a phosphocholine group with a negative charge"