"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine has a glycerol backbone, a phosphate group, two fatty acid chains, and an L-serine residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone linked to phosphate
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-[OX2]-[PX4](=[OX1])")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Check for the specific serine linkage pattern to the phosphate group, allowing for ionic forms.
    serine_phosphate_pattern = Chem.MolFromSmarts("[P]-[O]-[CH2X4]-[CHX4]([NX3])-[CX3](=[OX1])-[OX2]")
    serine_phosphate_pattern_ionic = Chem.MolFromSmarts("[P]([O-])-[O]-[CH2X4]-[CHX4]([NX3])-[CX3](=[OX1])-[OX2]")
    if not (mol.HasSubstructMatch(serine_phosphate_pattern) or mol.HasSubstructMatch(serine_phosphate_pattern_ionic)):
        return False, "No serine phosphate linkage found"
   

    # Check for two ester groups attached to glycerol at the correct position.
    # Updated patterns to explicitly require connection to the glycerol carbons.
    ester1_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4](-[OX2]-[CX3](=[OX1]))-[CH2X4]")  
    ester2_pattern = Chem.MolFromSmarts("[CH2X4](-[OX2]-[CX3](=[OX1]))-[CHX4]-[CH2X4]")
    ester_matches1 = mol.GetSubstructMatches(ester1_pattern)
    ester_matches2 = mol.GetSubstructMatches(ester2_pattern)
    
    if len(ester_matches1) + len(ester_matches2) != 2:
       return False, f"Found {len(ester_matches1) + len(ester_matches2)} esters attached to the glycerol, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # check the length of the fatty acid chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"


    # Count carbons, oxygens, nitrogen and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)


    if c_count < 10:
       return False, "Too few carbons"
    if o_count < 7:
       return False, "Too few oxygens"
    if p_count != 1:
       return False, "Must have exactly 1 phosphorus"
    if n_count != 1:
       return False, "Must have exactly 1 nitrogen"

    return True, "Contains glycerol backbone with phosphate group, serine residue, and 2 fatty acids"