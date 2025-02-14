"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:27230 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine is a glycerol backbone with a phosphoserine headgroup at the 3-position, and two fatty acid chains attached via ester bonds at the 1- and 2-positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glycerol backbone with phosphate at position 3 and correct stereochemistry
    # The pattern C[C@H](COP(...)) is crucial to determine the sn-3 configuration.
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](COP(=O)([OX1])[OX2])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 3-sn-glycerol-phosphate backbone found"

    # Check for the phosphoserine headgroup
    phosphoserine_pattern = Chem.MolFromSmarts("COP(=O)([OX1])OC[C@H]([NX3])C(=O)[OX2]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine headgroup found"
    
    # Check for two ester groups linked to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains - make sure they are at least medium length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"

    # Molecular weight should be somewhat high for a phospholipid
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
       return False, "Molecular weight too low for 3-sn-phosphatidyl-L-serine"

    # Count carbons, oxygens, nitrogen and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 15:
        return False, "Too few carbons for 3-sn-phosphatidyl-L-serine"
    if o_count < 8:
        return False, "Too few oxygens for 3-sn-phosphatidyl-L-serine"
    if n_count != 1:
        return False, "Must have exactly one nitrogen"
    if p_count != 1:
        return False, "Must have exactly one phosphorus"
    
    return True, "Contains 3-sn-glycerol backbone with two fatty acid chains and a phosphoserine headgroup"