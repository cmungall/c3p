"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:17855 1,2-diacyl-sn-glycero-3-phosphocholine

The conjugate base of a 1,2-diacyl-sn-glycero-3-phosphocholine compound formed by deprotonation of the phosphate OH group.
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

    # Check for phosphocholine group
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]  # Phosphorus atoms
    if not p_atoms:
        return False, "No phosphorus atom found"
    
    p_atom = p_atoms[0]
    o_atoms = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]  # Oxygen atoms
    if len(o_atoms) != 4:
        return False, "Incorrect number of oxygen atoms attached to phosphorus"
    
    negatively_charged_o = [o for o in o_atoms if o.GetFormalCharge() == -1]
    if not negatively_charged_o:
        return False, "No negatively charged oxygen atom found"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    return True, "Contains glycerol backbone with 2 fatty acid chains attached via ester bonds and a phosphocholine group with a negatively charged oxygen"